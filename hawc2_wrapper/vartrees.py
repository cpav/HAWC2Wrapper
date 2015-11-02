from numpy import zeros


class RotorVT(object):

    def __init__(self):

        self.hub_height = 0.0  # Hub height
        self.nblades = 3       # Number of blades
        self.tilt_angle = 0.0  # Tilt angle
        self.cone_angle = 0.0  # Cone angle
        self.diameter = 0.0    # Rotor diameter


class BladeVT(object):

    def __init__(self):
        self.length = 0.0   # blade length
        self.mass = 0.0     # blade mass
        self.I_x = 0.0      # first area moment of inertia
        self.I_y = 0.0      # Second area moment of inertia
        self.root_chord = 0.0  # Blade root chord
        self.max_chord = 0.0   # Blade maximum chord
        self.tip_chord = 0.0   # Blade tip chord
        self.airfoils = []     # List of airfoil names used on blade


class HubVT(object):

    def __init__(self):
        self.diameter = 0.0  # blade length
        self.mass = 0.0      # blade mass
        self.I_x = 0.0       # first area moment of inertia
        self.I_y = 0.0       # Second area moment of inertia
        self.CM = zeros(3)


class NacelleVT(object):

    def __init__(self):
        self.mass = 0.0  # blade mass
        self.I_x = 0.0   # first area moment of inertia
        self.I_y = 0.0   # Second area moment of inertia
        self.CM = zeros(3)
        self.diameter = 0.0


class ShaftVT(object):

    def __init__(self):
        self.mass = 0.0  # blade mass
        self.I_x = 0.0   # first area moment of inertia
        self.I_y = 0.0   # Second area moment of inertia
        self.CM = zeros(3)
        self.length = 0.0


class GeneratorVT(object):

    def __init__(self):
        self.mass = 0.0        # blade mass
        self.I_x = 0.0         # first area moment of inertia
        self.I_y = 0.0         # Second area moment of inertia
        self.power = 0.0       # Generator power
        self.max_torque = 0.0  # Maximum allowable generator torque
        self.efficiency = 0.0  # Generator efficiency


class TransmissionVT(object):

    def __init__(self):
        self.gear_ratio = 0.0  # Transmission gear ratio


class TowerVT(object):

    def __init__(self):
        self.height = 0.0           # Tower height
        self.bottom_diameter = 0.0  # Tower bottom diameter
        self.top_diameter = 0.0     # Tower bottom diameter
        self.mass = 0.0             # Tower mass


class BeamGeometryVT(object):

    def __init__(self):
        # Blade main axis accumulated curve length (n)
        self.s = zeros([1])
        self.x = zeros([1])
        self.y = zeros([1])
        self.z = zeros([1])
        self.rot_x = zeros([1])         # x-rotation angle (n)
        self.rot_y = zeros([1])         # y-rotation angle (n)
        self.rot_z = zeros([1])         # z-rotation angle (n)


class BladeGeometryVT(BeamGeometryVT):

    def __init__(self):
        super(BladeGeometryVT, self).__init__()
        self.chord = zeros([1])   # Blade chord (n)
        self.rthick = zeros([1])  # Blade relative thickness (n)
        self.athick = zeros([1])  # Blade absolute thickness (n)
        self.p_le = zeros([1])    # normalized distance along chord line from
                                  #  leading edge to main axis (n)


class DistributedLoadsVT(object):

    def __init__(self):
        self.s = zeros([1])   # locations for distributed loads
        self.Fn = zeros([1])  # force per unit length in normal direction
        self.Ft = zeros([1])  # force per unit length in tangential direction


class DistributedLoadsExtVT(DistributedLoadsVT):

    def __init__(self):
        super(DistributedLoadsVT, self).__init__()
        self.cn = zeros([1])  # Normal force coefficient along the blade
        self.ct = zeros([1])  # Tangential force coefficient along the blade
        self.cl = zeros([1])  # Lift force coefficient along the blade
        self.cd = zeros([1])  # Drag force coefficient along the blade
        self.cm = zeros([1])  # Moment force coefficient along the blade
        self.aoa = zeros([1])  # [deg] Angle of attack along the blade
        self.lfa = zeros([1])  # [deg] Local flow angle along the blade
        self.v_a = zeros([1])  # [m/s] axial velocity along the blade
        self.v_t = zeros([1])  # [m/s] tangential velocity along the blade
        self.v_r = zeros([1])  # [m/s] radial velocity along the blade
        self.lcp = zeros([1])  # Local power coefficient along the blade
        self.lct = zeros([1])  # Local power coefficient along the blade


class RotorLoadsVT(object):

    def __init__(self):
        self.T = 0.0  # [N] thrust
        self.Q = 0.0  # [N*m] torque
        self.P = 0.0  # [W] power

        self.CT = 0.0  # [N] thrust coefficient
        self.CQ = 0.0  # [N*m] torque coefficient
        self.CP = 0.0  # [W] power coefficient


class RotorLoadsArrayVT(object):

    def __init__(self):
        self.wsp = zeros([1])  # [m/s] Wind speeds
        self.rpm = zeros([1])  # [rpm] Rotor speed
        self.pitch = zeros([1])  # [deg] Pitch angle

        self.T = zeros([1])  # [N] thrust
        self.Q = zeros([1])  # [N*m] torque
        self.P = zeros([1])  # [W] power

        self.CT = zeros([1])  # thrust coefficient
        self.CQ = zeros([1])  # torque coefficient
        self.CP = zeros([1])  # power coefficient


class DistributedLoadsArrayVT(object):

    def __init__(self):
        self.loads_array = []


class BeamDisplacementsVT(object):

    def __init__(self):
        self.main_axis = zeros([1])
        self.rot_x = zeros([1])
        self.rot_y = zeros([1])
        self.rot_z = zeros([1])


class BeamDisplacementsArrayVT(object):

    def __init__(self):
        self.disps_array = []  # Array of blade displacements and rotations
        self.tip_pos = zeros([1])
        self.tip_rot = zeros([1])


class HubLoadsVT(object):

    def __init__(self):
        self.Fx = 0.0  # [N] x-force in wind-aligned coordinate system
        self.Fy = 0.0  # [N] y-force in wind-aligned coordinate system
        self.Fz = 0.0  # [N] z-force in wind-aligned coordinate system
        self.Mx = 0.0  # [N*m] x-moment in wind-aligned coordinate system
        self.My = 0.0  # [N*m] y-moment in wind-aligned coordinate system
        self.Mz = 0.0  # [N*m] z-moment in wind-aligned coordinate system


class HubLoadsArrayVT(object):

    def __init__(self):
        self.Fx = zeros([1])  # [N] x-force in wind-aligned coordinate system
        self.Fy = zeros([1])  # [N] y-force in wind-aligned coordinate system
        self.Fz = zeros([1])  # [N] z-force in wind-aligned coordinate system
        self.Mx = zeros([1])  # [Nm] x-moment in wind-aligned coordinate system
        self.My = zeros([1])  # [Nm] y-moment in wind-aligned coordinate system
        self.Mz = zeros([1])  # [Nm] z-moment in wind-aligned coordinate system


class RotorOperationalDataVT(object):

    def __init__(self):
        self.wsp = zeros([1])    # [m/s] Wind speed
        self.pitch = zeros([1])  # [deg] pitch angle
        self.rpm = zeros([1])    # rotational speed
