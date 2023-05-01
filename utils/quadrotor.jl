function dcm_from_mrp(p)
    p1,p2,p3 = p
    den = (p1^2 + p2^2 + p3^2 + 1)^2
    a = (4*p1^2 + 4*p2^2 + 4*p3^2 - 4)
    [
    (-((8*p2^2+8*p3^2)/den-1)*den)   (8*p1*p2 + p3*a)     (8*p1*p3 - p2*a);
    (8*p1*p2 - p3*a) (-((8*p1^2 + 8*p3^2)/den - 1)*den)   (8*p2*p3 + p1*a);
    (8*p1*p3 + p2*a)  (8*p2*p3 - p1*a)  (-((8*p1^2 + 8*p2^2)/den - 1)*den)
    ]/den
end
function skew(ω::Vector{T}) where {T}
    return [0 -ω[3] ω[2];
            ω[3] 0 -ω[1];
            -ω[2] ω[1] 0]
end

function dynamics_multi_drone(model::NamedTuple, x, u)
    # quadrotor dynamics with an MRP for attitude

    mass=model.mass
    J = model.J
    gravity= model.gravity
    L= model.L
    kf=model.kf
    km=model.km
    no_of_drones = model.no_of_drones

    
    result = []
    
    for drone = 1:no_of_drones
        u_single_drone = u[(drone-1)*4 + 1: (drone-1)*4 + 4]
        x_single_drone = x[(drone-1)*12 + 1: (drone-1)*12 + 12]

        r = x_single_drone[1:3]
        v = x_single_drone[4:6]
        p = x_single_drone[7:9]
        ω = x_single_drone[10:12]

        Q = dcm_from_mrp(p)
        
        w1 = u_single_drone[1]
        w2 = u_single_drone[2]
        w3 = u_single_drone[3]
        w4 = u_single_drone[4]

        F1 = max(0,kf*w1)
        F2 = max(0,kf*w2)
        F3 = max(0,kf*w3)
        F4 = max(0,kf*w4)
        F = [0., 0., F1+F2+F3+F4] #total rotor force in body frame

        M1 = km*w1
        M2 = km*w2
        M3 = km*w3
        M4 = km*w4
        τ = [L*(F2-F4), L*(F3-F1), (M1-M2+M3-M4)] #total rotor torque in body frame

        f = mass*gravity + Q*F # forces in world frame

        result = 
        [
            result
            v
            f/mass
            ((1+norm(p)^2)/4) *(   I + 2*(skew(p)^2 + skew(p))/(1+norm(p)^2)   )*ω
            J\(τ - cross(ω,J*ω))
        ]
    end
        
    return result
end

# function dynamics(model::NamedTuple,x,u)
    
#     for i = 1:no_of_drones
#         d = dynamics_single_drone(model::NamedTuple, x, u)
#     end
# #     @show d
#     @error ("dwec")
# end
function rk4(model,ode,x,u,dt)
    k1 = dt*ode(model,x, u)
    k2 = dt*ode(model,x + k1/2, u)
    k3 = dt*ode(model,x + k2/2, u)
    k4 = dt*ode(model,x + k3, u)
    x + (1/6)*(k1 + 2*k2 + 2*k3 + k4)
end
function vis_traj!(vis, name, X; R = 0.1, color = mc.RGBA(1.0, 0.0, 0.0, 1.0))
    for i = 1:(length(X)-1)
        a = X[i][1:3]
        b = X[i+1][1:3]
        cyl = mc.Cylinder(mc.Point(a...), mc.Point(b...), R)
        mc.setobject!(vis[name]["p"*string(i)], cyl, mc.MeshPhongMaterial(color=color))
    end
    for i = 1:length(X)
        a = X[i][1:3]
        sph = mc.HyperSphere(mc.Point(a...), R)
        mc.setobject!(vis[name]["s"*string(i)], sph, mc.MeshPhongMaterial(color=color))
    end
end

# Parameters:
# no_of_drones: integer
# Xsim: Vector of vectors. Each vector must be of size 12*no_of_drones
# dt: time step
function animate_quadrotor(no_of_drones, Xsim, dt)
    num_points = size(Xsim, 1)
    @assert num_points > 0
    
    num_states = size(Xsim[1], 1)
    @show num_states
    @assert num_states == no_of_drones * 12
    
    for k = 2:length(Xsim)
        num_states_curr = size(Xsim[k], 1)
        @assert num_states_curr == num_states
    end
    
    vis = mc.Visualizer()
    
    for drone_no = 1:no_of_drones                    
        robot_obj = mc.MeshFileGeometry(joinpath(@__DIR__,"quadrotor.obj"))
        drone_name = "drone" * string(drone_no)
        mc.setobject!(vis[Meta.parse(drone_name)], robot_obj)                 
    end

    anim = mc.Animation(floor(Int,1/dt))
    for drone_no = 1:no_of_drones
        for k = 1:length(Xsim)
            mc.atframe(anim, k) do
                                    
                r = Xsim[k][((drone_no - 1) * 12) + 1 : ((drone_no - 1) * 12) + 3]
                p = Xsim[k][((drone_no - 1) * 12) + 7 : ((drone_no - 1) * 12) + 9]
                drone_name = "drone" * string(drone_no)
                mc.settransform!(vis[Meta.parse(drone_name)], mc.compose(mc.Translation(r),mc.LinearMap(1.5*(dcm_from_mrp(p)))))
            end
        end
            mc.setanimation!(vis, anim)
    end
    return (mc.render(vis))
end