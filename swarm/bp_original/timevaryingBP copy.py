import torch
import numpy as np
import time
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

# Check for GPU
DEVICE = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
print(f"Using device: {DEVICE}")

class BPSDF_2D_Model_With_Grad:
    def __init__(self, n_func, domain_min, domain_max):
        self.n_func = n_func
        self.domain_min = torch.tensor(domain_min, dtype=torch.float32, device=DEVICE)
        self.domain_max = torch.tensor(domain_max, dtype=torch.float32, device=DEVICE)
        self.device = DEVICE
        self.comb_matrix = self._precompute_comb()

    def _precompute_comb(self):
        n = self.n_func - 1
        k = torch.arange(self.n_func, device=self.device, dtype=torch.float32)
        return torch.exp(torch.lgamma(torch.tensor(n + 1, device=self.device)) - 
                         torch.lgamma(k + 1) - 
                         torch.lgamma(torch.tensor(n + 1, device=self.device) - k))

    def build_bernstein_t(self, t, use_derivative=False):
        t = torch.clamp(t, min=1e-6, max=1 - 1e-6)
        n = self.n_func - 1
        i = torch.arange(self.n_func, device=self.device)
        comb = self.comb_matrix
        if t.ndim == 0: t = t.unsqueeze(0)
        term1 = (1 - t).unsqueeze(-1) ** (n - i)
        term2 = t.unsqueeze(-1) ** i
        phi = comb * term1 * term2

        if not use_derivative:
            return phi.float(), None
        else:
            term1_d = i * t.unsqueeze(-1)**torch.clamp(i-1, min=0) * (1-t).unsqueeze(-1)**(n-i)
            term2_d = (n-i) * t.unsqueeze(-1)**i * (1-t).unsqueeze(-1)**torch.clamp(n-i-1, min=0)
            dphi = comb * (term1_d - term2_d)
            return phi.float(), dphi.float()

    def build_basis_function_from_points_2d(self, p):
        p = torch.atleast_2d(p).float().to(self.device)
        p_normalized = ((p - self.domain_min) / (self.domain_max - self.domain_min))
        phi_x, _ = self.build_bernstein_t(p_normalized[:, 0], use_derivative=False)
        phi_y, _ = self.build_bernstein_t(p_normalized[:, 1], use_derivative=False)
        phi_xy = torch.einsum("bi,bj->bij", phi_x, phi_y).reshape(p.shape[0], self.n_func**2)
        return phi_xy

    def build_gradient_basis_from_points_2d(self, p):
        p = torch.atleast_2d(p).float().to(self.device)
        domain_span = self.domain_max - self.domain_min
        p_normalized = ((p - self.domain_min) / (domain_span))
        phi_x, dphi_x = self.build_bernstein_t(p_normalized[:, 0], use_derivative=True)
        phi_y, dphi_y = self.build_bernstein_t(p_normalized[:, 1], use_derivative=True)
        grad_x_basis = torch.einsum("bi,bj->bij", dphi_x, phi_y).reshape(p.shape[0], self.n_func**2) / domain_span[0]
        grad_y_basis = torch.einsum("bi,bj->bij", phi_x, dphi_y).reshape(p.shape[0], self.n_func**2) / domain_span[1]
        grad_basis_stacked = torch.cat([grad_x_basis, grad_y_basis], dim=0)
        return grad_basis_stacked
        
    def predict(self, points, weights):
        # Returns a pure PyTorch Tensor of shape (N,)
        with torch.no_grad():
            psi = self.build_basis_function_from_points_2d(points)
            return (psi @ weights.view(-1, 1)).squeeze()

    def gradient(self, points, weights):
        # Returns a pure PyTorch Tensor of shape (N, 2)
        with torch.no_grad():
            weights_flat = weights.view(-1, 1)
            grad_psi = self.build_gradient_basis_from_points_2d(points)
            num_points = points.shape[0]
            grad_x = grad_psi[:num_points] @ weights_flat
            grad_y = grad_psi[num_points:] @ weights_flat
            return torch.cat([grad_x, grad_y], dim=1)

# --- Vectorized Helper Functions ---
def get_scene_sdf_vectorized(points, centers, radius):
    diff = points[:, np.newaxis, :] - centers[np.newaxis, :, :] 
    dists = np.linalg.norm(diff, axis=2) 
    min_dists = np.min(dists, axis=1) 
    return min_dists - radius

def generate_training_data_batch(num_surface_points, num_random_points, centers_tensor, radius, domain_min, domain_max):
    SAFETY_MARGIN = 0.0  
    centers = centers_tensor.cpu().numpy()
    num_circles = centers.shape[0]
    
    if num_circles == 0:
        return torch.empty(0, 2), torch.empty(0), torch.empty(0, 2), torch.empty(0, 2)

    # Surface points
    points_per_circle = max(1, num_surface_points // num_circles)
    angles = np.linspace(0, 2 * np.pi, points_per_circle, endpoint=False)
    cos_sin = np.vstack([np.cos(angles), np.sin(angles)]).T * radius
    
    grad_points = (centers[:, np.newaxis, :] + cos_sin[np.newaxis, :, :]).reshape(-1, 2)
    grad_normals = np.tile(cos_sin / radius, (num_circles, 1))
    
    surface_sdfs = np.zeros(len(grad_points)) - SAFETY_MARGIN

    # Random points
    random_points = np.random.uniform(low=domain_min, high=domain_max, size=(num_random_points, 2))
    random_sdfs = get_scene_sdf_vectorized(random_points, centers, radius) - SAFETY_MARGIN

    all_points = np.vstack([grad_points, random_points])
    all_sdfs = np.concatenate([surface_sdfs, random_sdfs])
    
    return (torch.tensor(all_points, dtype=torch.float32, device=DEVICE),
            torch.tensor(all_sdfs, dtype=torch.float32, device=DEVICE),
            torch.tensor(grad_points, dtype=torch.float32, device=DEVICE),
            torch.tensor(grad_normals, dtype=torch.float32, device=DEVICE))

def train_batch_with_gradient(model, points, sdfs, grad_points, grad_normals):
    psi_sdf = model.build_basis_function_from_points_2d(points)
    psi_grad = model.build_gradient_basis_from_points_2d(grad_points)
    grad_weight = 0.5
    psi_grad *= grad_weight
    psi_all = torch.cat([psi_sdf, psi_grad], dim=0)
    target_f = torch.cat([sdfs, grad_weight * grad_normals.T.flatten()], dim=0)
    
    lambda_reg = 1e-4
    psi_t = psi_all.T
    A = psi_t @ psi_all
    A.diagonal().add_(lambda_reg)
    b = psi_t @ target_f
    return torch.linalg.solve(A, b)

# --- BATCH KKT SOLVER (Pure PyTorch) ---
def solve_cbf_analytic_batch(u_des, A, b, max_speed):
    """
    Solves A * u >= b for an entire batch of agents simultaneously.
    u_des: (N, 2), A: (N, 2), b: (N,)
    """
    constraint_val = (A * u_des).sum(dim=1) - b
    norm_A_sq = (A * A).sum(dim=1)
    
    # Avoid division by zero
    safe_norm_A_sq = torch.where(norm_A_sq < 1e-8, torch.ones_like(norm_A_sq), norm_A_sq)
    
    lam = (b - (A * u_des).sum(dim=1)) / safe_norm_A_sq
    lam = torch.clamp(lam, min=0.0) # Apply only if constraint is violated
    
    u_safe = u_des + lam.unsqueeze(1) * A
    
    # Clip to max speed
    speeds = torch.norm(u_safe, dim=1)
    clip_mask = speeds > max_speed
    if clip_mask.any():
        u_safe[clip_mask] = (u_safe[clip_mask] / speeds[clip_mask].unsqueeze(1)) * max_speed
        
    return u_safe

# --- FULL VECTORIZED SWARM CONTROL ---
def get_safe_control_batch(pos, vels, goals, radius, other_model, weights_other, weights_dot_other, gamma, max_speed):
    N = pos.shape[0]
    
    # 1. Nominal Desired Velocity (Vectorized)
    diff = goals - pos
    dist_goal = torch.norm(diff, dim=1, keepdim=True)
    u_des = torch.where(dist_goal > 0.1, (diff / dist_goal) * max_speed, torch.zeros_like(diff))
        
    # --- 2. Inter-Swarm Constraint (Level Set) ---
    sdf_val = other_model.predict(pos, weights_other) # Shape: (N,)
    min_h = sdf_val - radius
    
    grad_h = other_model.gradient(pos, weights_other) # Shape: (N, 2)
    psi_crit = other_model.build_basis_function_from_points_2d(pos)
    dh_dt_shape = (psi_crit @ weights_dot_other.view(-1, 1)).squeeze() # Shape: (N,)
    
    A_inter = grad_h
    b_inter = -gamma * min_h - dh_dt_shape
    
    # --- 3. Intra-Swarm Constraints (Avoid Teammates) ---
    dist_matrix = torch.cdist(pos, pos)
    dist_matrix.fill_diagonal_(float('inf'))
    
    min_intra_dist, min_intra_idx = torch.min(dist_matrix, dim=1)
    h_intra = min_intra_dist - (2 * radius)
    
    pos_teammate = pos[min_intra_idx]
    vel_teammate = vels[min_intra_idx]
    
    grad_intra = (pos - pos_teammate) / min_intra_dist.unsqueeze(1)
    grad_intra = torch.where(min_intra_dist.unsqueeze(1) > 1e-6, grad_intra, torch.tensor([1.0, 0.0], device=DEVICE))
    b_intra = -gamma * h_intra + (grad_intra * vel_teammate).sum(dim=1)
    
    # --- 4. Most Binding Constraint Selection ---
    # Select whichever constraint (inter vs intra) has the smallest h
    inter_is_binding = min_h < h_intra
    
    best_A = torch.where(inter_is_binding.unsqueeze(1), A_inter, grad_intra)
    best_b = torch.where(inter_is_binding, b_inter, b_intra)
        
    # --- 5. Batch KKT Projection ---
    u_safe = solve_cbf_analytic_batch(u_des, best_A, best_b, max_speed)
    return u_safe

if __name__ == '__main__':
    np.random.seed(42)
    torch.manual_seed(42)

    N_func = 24
    domain_size = 4.5
    domain_min = np.array([-domain_size, -domain_size])
    domain_max = np.array([domain_size, domain_size])
    
    model_A = BPSDF_2D_Model_With_Grad(n_func=N_func, domain_min=domain_min, domain_max=domain_max)
    model_B = BPSDF_2D_Model_With_Grad(n_func=N_func, domain_min=domain_min, domain_max=domain_max)
    
    # Swarm A (Blue): Tight V-formation, moves Top -> Bottom
    pos_A = torch.tensor([[0.0, 3.0], [-0.8, 3.8], [0.8, 3.8]], dtype=torch.float32, device=DEVICE)
    goals_A = torch.tensor([[0.0, -3.5], [-0.8, -3.5], [0.8, -3.5]], dtype=torch.float32, device=DEVICE)
    vels_A = torch.zeros_like(pos_A)
    radius_A = 0.4
    
    # Swarm B (Red): Tight V-formation, moves Bottom -> Top (with a 0.1 x-offset)
    pos_B = torch.tensor([[0.1, -3.0], [-0.7, -3.8], [0.9, -3.8]], dtype=torch.float32, device=DEVICE)
    goals_B = torch.tensor([[0.1, 3.5], [-0.7, 3.5], [0.9, 3.5]], dtype=torch.float32, device=DEVICE)
    vels_B = torch.zeros_like(pos_B)
    radius_B = 0.4

    num_frames = 1000
    dt = 0.02 
    gamma = 40.0 
    max_speed = 2.0

    plt.ion()
    fig, ax1 = plt.subplots(1, 1, figsize=(8, 8))
    fig.canvas.manager.set_window_title("GPU Vectorized Level-Set CBF Intersection")
    fig.tight_layout()

    vis_res = 120
    x_vis, y_vis = np.meshgrid(np.linspace(domain_min[0], domain_max[0], vis_res),
                              np.linspace(domain_min[1], domain_max[1], vis_res))
    points_to_visualize = torch.tensor(np.vstack([x_vis.ravel(), y_vis.ravel()]).T, dtype=torch.float32, device=DEVICE)

    # --- EMA Setup ---
    weights_A_prev, weights_A_dot = None, None
    weights_B_prev, weights_B_dot = None, None
    alpha_ema = 0.2 

    print(f"--- Starting GPU Accelerated Dual Swarm Simulation on {DEVICE} ---")

    for frame in range(num_frames):
        t_start = time.perf_counter()

        # 1) Generate Data & Train
        pts_A, sdf_A, grad_pts_A, grad_n_A = generate_training_data_batch(6000, 15000, pos_A, radius_A, -domain_size, domain_size)
        weights_A = train_batch_with_gradient(model_A, pts_A, sdf_A, grad_pts_A, grad_n_A)

        pts_B, sdf_B, grad_pts_B, grad_n_B = generate_training_data_batch(6000, 15000, pos_B, radius_B, -domain_size, domain_size)
        weights_B = train_batch_with_gradient(model_B, pts_B, sdf_B, grad_pts_B, grad_n_B)

        # 2) Track Shape Deformation
        if weights_A_prev is None:
            weights_A_prev, weights_B_prev = weights_A.clone(), weights_B.clone()
            weights_A_dot, weights_B_dot = torch.zeros_like(weights_A), torch.zeros_like(weights_B)
        else:
            weights_A_dot = (alpha_ema * ((weights_A - weights_A_prev) / dt)) + ((1.0 - alpha_ema) * weights_A_dot)
            weights_B_dot = (alpha_ema * ((weights_B - weights_B_prev) / dt)) + ((1.0 - alpha_ema) * weights_B_dot)
            weights_A_prev, weights_B_prev = weights_A.clone(), weights_B.clone()

        # 3) Vectorized GPU Control Output
        vels_A_new = get_safe_control_batch(pos_A, vels_A, goals_A, radius_A, model_B, weights_B, weights_B_dot, gamma, max_speed)
        vels_B_new = get_safe_control_batch(pos_B, vels_B, goals_B, radius_B, model_A, weights_A, weights_A_dot, gamma, max_speed)

        # 4) Apply Controls
        vels_A = vels_A_new
        pos_A += vels_A * dt
        
        vels_B = vels_B_new
        pos_B += vels_B * dt

        # --- Plotting (Only move to CPU here) ---
        ax1.clear()
        
        z_A = model_A.predict(points_to_visualize, weights_A).cpu().numpy().reshape(vis_res, vis_res)
        z_B = model_B.predict(points_to_visualize, weights_B).cpu().numpy().reshape(vis_res, vis_res)
        
        ax1.contourf(x_vis, y_vis, z_A, levels=[-100, 0], colors=['royalblue'], alpha=0.3)
        ax1.contour(x_vis, y_vis, z_A, levels=[0], colors='blue', linewidths=2)
        
        ax1.contourf(x_vis, y_vis, z_B, levels=[-100, 0], colors=['crimson'], alpha=0.3)
        ax1.contour(x_vis, y_vis, z_B, levels=[0], colors='red', linewidths=2)

        pos_A_np = pos_A.cpu().numpy()
        goals_A_np = goals_A.cpu().numpy()
        for i in range(pos_A_np.shape[0]):
            ax1.add_patch(Circle(pos_A_np[i], radius_A, color='midnightblue', fill=True, alpha=0.9))
            ax1.plot(goals_A_np[i, 0], goals_A_np[i, 1], 'x', color='blue', ms=8)

        pos_B_np = pos_B.cpu().numpy()
        goals_B_np = goals_B.cpu().numpy()
        for i in range(pos_B_np.shape[0]):
            ax1.add_patch(Circle(pos_B_np[i], radius_B, color='darkred', fill=True, alpha=0.9))
            ax1.plot(goals_B_np[i, 0], goals_B_np[i, 1], 'x', color='red', ms=8)

        ax1.set_title(f'GPU Vectorized KKT - Frame {frame+1}')
        ax1.set_xlim(domain_min[0], domain_max[0]); ax1.set_ylim(domain_min[1], domain_max[1])
        ax1.set_aspect('equal')
        
        plt.pause(0.001)
        t_end = time.perf_counter()
        print(f"Frame {frame+1:03d} | Compute: {(t_end - t_start):.4f}s")

    plt.ioff()
    plt.show()