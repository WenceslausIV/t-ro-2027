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
        with torch.no_grad():
            psi = self.build_basis_function_from_points_2d(points)
            return (psi @ weights.view(-1, 1)).squeeze()

    def gradient(self, points, weights):
        with torch.no_grad():
            weights_flat = weights.view(-1, 1)
            grad_psi = self.build_gradient_basis_from_points_2d(points)
            num_points = points.shape[0]
            grad_x = grad_psi[:num_points] @ weights_flat
            grad_y = grad_psi[num_points:] @ weights_flat
            return torch.cat([grad_x, grad_y], dim=1)

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

    points_per_circle = max(1, num_surface_points // num_circles)
    angles = np.linspace(0, 2 * np.pi, points_per_circle, endpoint=False)
    cos_sin = np.vstack([np.cos(angles), np.sin(angles)]).T * radius
    
    grad_points_unfiltered = (centers[:, np.newaxis, :] + cos_sin[np.newaxis, :, :]).reshape(-1, 2)
    grad_normals_unfiltered = np.tile(cos_sin / radius, (num_circles, 1))
    
    true_sdfs = get_scene_sdf_vectorized(grad_points_unfiltered, centers, radius)
    outer_mask = true_sdfs >= -1e-4 
    
    grad_points = grad_points_unfiltered[outer_mask]
    grad_normals = grad_normals_unfiltered[outer_mask]
    surface_sdfs = np.zeros(len(grad_points)) - SAFETY_MARGIN

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

def solve_cbf_analytic_batch(u_des, A, b, max_speed):
    constraint_val = (A * u_des).sum(dim=1) - b
    norm_A_sq = (A * A).sum(dim=1)
    
    safe_norm_A_sq = torch.where(norm_A_sq < 1e-8, torch.ones_like(norm_A_sq), norm_A_sq)
    
    lam = (b - (A * u_des).sum(dim=1)) / safe_norm_A_sq
    lam = torch.clamp(lam, min=0.0) 
    
    u_safe = u_des + lam.unsqueeze(1) * A
    
    speeds = torch.norm(u_safe, dim=1)
    clip_mask = speeds > max_speed
    if clip_mask.any():
        u_safe[clip_mask] = (u_safe[clip_mask] / speeds[clip_mask].unsqueeze(1)) * max_speed
        
    return u_safe

def apply_intra_swarm_cbf(u_des, pos, radius, gamma_inner, max_speed, alpha_lse=20.0):
    """
    Basic circle-to-circle CBF for intra-swarm (same-color) collision avoidance.

    Barrier function for robot pair (i, j):
        h_ij = ||p_i - p_j||^2 - (2r)^2   (>= 0 when separated, < 0 when overlapping)

    CBF constraint for robot i (treating neighbor j's velocity as unknown → 0):
        2*(p_i - p_j)^T * u_i  >=  -gamma_inner * h_ij

    Multiple neighbors are aggregated via LSE soft-minimum (same as inter-swarm CBF),
    yielding a single linear constraint per robot that is solved analytically.
    """
    N = pos.shape[0]
    if N <= 1:
        return u_des

    # Pairwise difference vectors and squared distances
    diff = pos.unsqueeze(1) - pos.unsqueeze(0)          # (N, N, 2)  diff[i,j] = p_i - p_j
    dist_sq = torch.sum(diff ** 2, dim=2)               # (N, N)
    h_mat = dist_sq - (2.0 * radius) ** 2              # (N, N)  h_ij

    # Strip diagonal (self-pairs) → (N, N-1) views
    eye_mask = ~torch.eye(N, dtype=torch.bool, device=pos.device)
    h_neighbors    = h_mat[eye_mask].view(N, N - 1)         # (N, N-1)
    diff_neighbors = diff[eye_mask].view(N, N - 1, 2)       # (N, N-1, 2)

    # LSE smooth minimum over all (N-1) neighbors
    min_h_smooth = -(1.0 / alpha_lse) * torch.logsumexp(-alpha_lse * h_neighbors, dim=1)  # (N,)
    weights_lse  = torch.softmax(-alpha_lse * h_neighbors, dim=1)                          # (N, N-1)

    # Aggregated constraint direction:  A_i = sum_j w_j * 2*(p_i - p_j)
    A_inner = 2.0 * torch.sum(weights_lse.unsqueeze(2) * diff_neighbors, dim=1)  # (N, 2)
    b_inner = -gamma_inner * min_h_smooth                                          # (N,)

    return solve_cbf_analytic_batch(u_des, A_inner, b_inner, max_speed)


def get_adaptive_level_set(model, weights, pos, radius, n_radii=6, n_angles=32, r_max=1.1):
    """
    Compute the adaptive level set that fully encloses all robot disks.

    Key insight: sampling only UP TO r = radius (true boundary) means every sampled
    point has true SDF <= 0.  If the BP model underestimates the boundary (assigns
    negative phi to boundary points due to fitting error), then
        l = max(phi over disk) < 0
    and the isosurface phi = l falls INSIDE the robot — enclosure fails.

    Fix: extend sampling to r_max = 1.1 * radius (slightly outside each robot).
    Points outside a robot have true SDF > 0, so the BP model gives phi > 0 there
    even with approximation errors.  This forces l > 0, placing the isosurface
    strictly outside all robots and guaranteeing full enclosure.
    """
    # Build angle samples — skip duplicate endpoint
    angles = torch.linspace(0, 2 * np.pi, n_angles + 1, device=DEVICE)[:-1]
    cos_sin = torch.stack([torch.cos(angles), torch.sin(angles)], dim=1)  # (n_angles, 2)

    all_pts = [pos]  # robot centres (r = 0, true SDF = -radius)
    # Rings from interior through boundary and slightly beyond (r_max > 1.0)
    # The exterior ring (r > 1.0) is what forces l > 0 and guarantees enclosure.
    for r_scale in torch.linspace(0.0, r_max, n_radii, device=DEVICE)[1:]:
        ring_pts = pos.unsqueeze(1) + (r_scale * radius) * cos_sin.unsqueeze(0)  # (N, n_angles, 2)
        all_pts.append(ring_pts.view(-1, 2))

    all_pts_flat = torch.cat(all_pts, dim=0)
    vals = model.predict(all_pts_flat, weights)
    return torch.max(vals)

def get_safe_control_exact_adaptive(pos, goals, radius,
                                    my_model, my_weights, l_my, l_my_dot,
                                    other_model, weights_other, weights_dot_other, l_other, l_dot_other,
                                    gamma, max_speed, alpha_lse=20.0):
    """
    CBF-QP safe control with fully adaptive level sets.

    h = phi_other(p) - l_other   where p = projection of robot surface onto MY level set phi_my = l_my

    Total dh/dt (autonomous, non-control part) has three contributions:
      1. dh/dt_weights : phi_other changes because weights_other change  → psi_other · w_dot_other
      2. dh/dt_l_other : l_other changes explicitly                      → -l_dot_other
      3. dh/dt_l_my    : l_my changes, moving the projected points p
                         dp/dt = l_my_dot * ∇phi_my / |∇phi_my|²
                         → ∇phi_other(p) · dp/dt
    All three are folded into b_cbf so the QP remains a simple linear programme.
    """
    N = pos.shape[0]

    diff = goals - pos
    dist_goal = torch.norm(diff, dim=1, keepdim=True)
    u_des = torch.where(dist_goal > 0.1, (diff / dist_goal) * max_speed, torch.zeros_like(diff))

    angles = torch.linspace(0, 2 * np.pi, 16, device=DEVICE)
    cos_sin = torch.stack([torch.cos(angles), torch.sin(angles)], dim=1)
    pts = pos.unsqueeze(1) + radius * cos_sin.unsqueeze(0)
    pts_flat = pts.view(-1, 2)  # (N*16, 2)

    # Project pts onto MY level-set isosurface phi_my = l_my
    my_sdf_val = my_model.predict(pts_flat, my_weights) - l_my   # (N*16,)
    my_grad = my_model.gradient(pts_flat, my_weights)             # (N*16, 2)
    my_grad_norm_sq = torch.sum(my_grad ** 2, dim=1)              # (N*16,)
    safe_grad_norm_sq = torch.where(my_grad_norm_sq < 1e-6,
                                    torch.ones_like(my_grad_norm_sq),
                                    my_grad_norm_sq)

    pts_exact_flat = pts_flat - (my_sdf_val / safe_grad_norm_sq).unsqueeze(1) * my_grad  # (N*16, 2)

    # Evaluate h = phi_other(p) - l_other at projected points
    h_vals_flat = other_model.predict(pts_exact_flat, weights_other) - l_other  # (N*16,)
    h_vals = h_vals_flat.view(N, 16)

    # --- LSE smooth minimum ---
    min_h_smooth = -(1.0 / alpha_lse) * torch.logsumexp(-alpha_lse * h_vals, dim=1)  # (N,)
    weights_lse = torch.softmax(-alpha_lse * h_vals, dim=1)                           # (N, 16)

    # Gradients and time derivatives of h at projected points
    grad_h_all = other_model.gradient(pts_exact_flat, weights_other).view(N, 16, 2)  # (N, 16, 2)
    psi_crit_all = other_model.build_basis_function_from_points_2d(pts_exact_flat)
    dh_dt_weights = (psi_crit_all @ weights_dot_other.view(-1, 1)).view(N, 16)       # (N, 16)  term①

    # Contribution from MY level-set speed (term③):
    #   dp/dt_l_my = l_my_dot * ∇phi_my / |∇phi_my|²
    #   dh/dt_l_my = ∇phi_other(p) · dp/dt_l_my
    dp_coeff = l_my_dot / safe_grad_norm_sq                    # (N*16,) scalar per projected point
    dh_dt_l_my_flat = torch.sum(
        grad_h_all.view(N * 16, 2) * (dp_coeff.unsqueeze(1) * my_grad), dim=1
    )                                                          # (N*16,)
    dh_dt_l_my = dh_dt_l_my_flat.view(N, 16)                  # (N, 16)  term③

    # LSE-weighted smooth aggregation
    grad_h_smooth = torch.sum(weights_lse.unsqueeze(2) * grad_h_all, dim=1)         # (N, 2)
    dh_dt_smooth  = torch.sum(weights_lse * (dh_dt_weights + dh_dt_l_my), dim=1)   # (N,)  ①+③

    # CBF constraint:  A·u >= b
    #   where b = -γ·h - (dh/dt_weights + dh/dt_l_my) + l_dot_other   (term② is +l_dot_other)
    A_cbf = grad_h_smooth
    b_cbf = -gamma * min_h_smooth - (dh_dt_smooth - l_dot_other)

    u_safe = solve_cbf_analytic_batch(u_des, A_cbf, b_cbf, max_speed)
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
    
    pos_A = torch.tensor([
        [-0.5, 3.0], [-1.3, 3.0], [0.3, 3.0], 
        [-0.9, 3.8], [-0.1, 3.8]              
    ], dtype=torch.float32, device=DEVICE)
    goals_A = torch.tensor([
        [-0.5, -3.5], [-1.3, -3.5], [0.3, -3.5],
        [-0.9, -4.3], [-0.1, -4.3]
    ], dtype=torch.float32, device=DEVICE)
    radius_A = 0.4
    
    pos_B = torch.tensor([
        [0.5, -3.0], [-0.3, -3.0], [1.3, -3.0], 
        [0.1, -3.8], [0.9, -3.8]               
    ], dtype=torch.float32, device=DEVICE)
    goals_B = torch.tensor([
        [0.5, 3.5], [-0.3, 3.5], [1.3, 3.5],
        [0.1, 4.3], [0.9, 4.3]
    ], dtype=torch.float32, device=DEVICE)
    radius_B = 0.4

    num_frames = 1000
    dt = 0.02
    gamma       = 60.0   # inter-swarm BP-CBF gain
    gamma_inner = 30.0   # intra-swarm circle-CBF gain (same-color robots)
    max_speed   = 2.0

    plt.ion()
    fig, ax1 = plt.subplots(1, 1, figsize=(8, 8))
    fig.canvas.manager.set_window_title("BP-CBF Inter-Swarm + Circle-CBF Intra-Swarm")
    fig.tight_layout()

    vis_res = 120
    x_vis, y_vis = np.meshgrid(np.linspace(domain_min[0], domain_max[0], vis_res),
                              np.linspace(domain_min[1], domain_max[1], vis_res))
    points_to_visualize = torch.tensor(np.vstack([x_vis.ravel(), y_vis.ravel()]).T, dtype=torch.float32, device=DEVICE)

    weights_A_prev, weights_A_dot = None, None
    weights_B_prev, weights_B_dot = None, None
    l_A_prev, l_B_prev = None, None
    l_A_dot, l_B_dot = 0.0, 0.0
    alpha_ema = 0.2 

    print(f"--- Starting Outer-Hull Filtered Simulation on {DEVICE} ---")

    for frame in range(num_frames):
        t_start = time.perf_counter()

        pts_A, sdf_A, grad_pts_A, grad_n_A = generate_training_data_batch(6000, 15000, pos_A, radius_A, -domain_size, domain_size)
        weights_A = train_batch_with_gradient(model_A, pts_A, sdf_A, grad_pts_A, grad_n_A)

        pts_B, sdf_B, grad_pts_B, grad_n_B = generate_training_data_batch(6000, 15000, pos_B, radius_B, -domain_size, domain_size)
        weights_B = train_batch_with_gradient(model_B, pts_B, sdf_B, grad_pts_B, grad_n_B)

        l_A = get_adaptive_level_set(model_A, weights_A, pos_A, radius_A)
        l_B = get_adaptive_level_set(model_B, weights_B, pos_B, radius_B)

        if weights_A_prev is None:
            weights_A_prev, weights_B_prev = weights_A.clone(), weights_B.clone()
            weights_A_dot, weights_B_dot = torch.zeros_like(weights_A), torch.zeros_like(weights_B)
            l_A_prev, l_B_prev = l_A.clone(), l_B.clone()
        else:
            weights_A_dot = (alpha_ema * ((weights_A - weights_A_prev) / dt)) + ((1.0 - alpha_ema) * weights_A_dot)
            weights_B_dot = (alpha_ema * ((weights_B - weights_B_prev) / dt)) + ((1.0 - alpha_ema) * weights_B_dot)
            
            l_A_dot = (alpha_ema * ((l_A - l_A_prev) / dt)) + ((1.0 - alpha_ema) * l_A_dot)
            l_B_dot = (alpha_ema * ((l_B - l_B_prev) / dt)) + ((1.0 - alpha_ema) * l_B_dot)

            weights_A_prev, weights_B_prev = weights_A.clone(), weights_B.clone()
            l_A_prev, l_B_prev = l_A.clone(), l_B.clone()

        # --- Stage 1: inter-swarm collision avoidance (BP-SDF CBF) ---
        vels_A_inter = get_safe_control_exact_adaptive(
            pos_A, goals_A, radius_A,
            model_A, weights_A, l_A, l_A_dot,
            model_B, weights_B, weights_B_dot, l_B, l_B_dot,
            gamma, max_speed
        )
        vels_B_inter = get_safe_control_exact_adaptive(
            pos_B, goals_B, radius_B,
            model_B, weights_B, l_B, l_B_dot,
            model_A, weights_A, weights_A_dot, l_A, l_A_dot,
            gamma, max_speed
        )

        # --- Stage 2: intra-swarm collision avoidance (basic circle CBF, same color) ---
        vels_A_new = apply_intra_swarm_cbf(vels_A_inter, pos_A, radius_A, gamma_inner, max_speed)
        vels_B_new = apply_intra_swarm_cbf(vels_B_inter, pos_B, radius_B, gamma_inner, max_speed)

        pos_A += vels_A_new * dt
        pos_B += vels_B_new * dt

        ax1.clear()
        
        z_A = model_A.predict(points_to_visualize, weights_A).cpu().numpy().reshape(vis_res, vis_res)
        z_B = model_B.predict(points_to_visualize, weights_B).cpu().numpy().reshape(vis_res, vis_res)
        
        l_A_np = l_A.item()
        l_B_np = l_B.item()

        ax1.contourf(x_vis, y_vis, z_A, levels=[-100, l_A_np], colors=['royalblue'], alpha=0.3)
        ax1.contour(x_vis, y_vis, z_A, levels=[l_A_np], colors='blue', linewidths=2)
        
        ax1.contourf(x_vis, y_vis, z_B, levels=[-100, l_B_np], colors=['crimson'], alpha=0.3)
        ax1.contour(x_vis, y_vis, z_B, levels=[l_B_np], colors='red', linewidths=2)

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

        ax1.set_title(f'Outer-Hull Filtered Level-Set Interaction (l_A:{l_A_np:.3f}, l_B:{l_B_np:.3f})')
        ax1.set_xlim(domain_min[0], domain_max[0]); ax1.set_ylim(domain_min[1], domain_max[1])
        ax1.set_aspect('equal')
        
        plt.pause(0.001)
        t_end = time.perf_counter()
        print(f"Frame {frame+1:03d} | l_A: {l_A_np:.4f} | l_B: {l_B_np:.4f} | Compute: {(t_end - t_start):.4f}s")

    plt.ioff()
    plt.show()