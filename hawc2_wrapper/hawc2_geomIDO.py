
import numpy as np
from scipy.interpolate import pchip

from hawc2_vartrees import HAWC2BladeGeometry, HAWC2BeamStructure
from vartrees import BladeGeometryVT

from PGL.main.distfunc import distfunc


class HAWC2GeometryBuilder(object):
    """
    parameters
    ----------
    bladegeom: BladeGeometryVT

    c12axis_init: array

    interp_from_htc: bool
        Interpolate blade onto the distribution defined in the htc master file

    blade_ni_span: int
        spanwise distribution of blade planform

    blade_length: float

    hub_radius: float

    returns
    -------
    blade_ae: HAWC2BladeGeometry

    c12axis: array

    """
    def __init__(self, **kwargs):
        super(HAWC2GeometryBuilder, self).__init__()

        self.bladegeom = BladeGeometryVT()
        self.blade_ae = HAWC2BladeGeometry()

        for k, w in kwargs.iteritems():
            try:
                setattr(self, k, w)
            except:
                pass

    def execute(self):

        if self.interp_from_htc:
            c12axis = self.c12axis_init
        else:
            c12axis = self.calculate_c12axis()

        ds_root = 1. / self.blade_ni_span
        ds_tip = 1. / self.blade_ni_span / 3.
        dist = np.array([[0., ds_root, 1], [1., ds_tip, self.blade_ni_span]])

        x = distfunc(dist)
        self.c12axis = np.zeros((x.shape[0], 4))
        for i in range(4):
            tck = pchip(c12axis[:, 2], c12axis[:, i])
            self.c12axis[:, i] = tck(x)

        # scale main axis according to radius
        self.c12axis[:, :3] *= self.blade_length

        self.blade_ae.radius = self.blade_length + self.hub_radius
        l = ((self.c12axis[1:, 0]-self.c12axis[:-1, 0])**2 + 
             (self.c12axis[1:, 1]-self.c12axis[:-1, 1])**2 +
             (self.c12axis[1:, 2]-self.c12axis[:-1, 2])**2)**.5
        self.blade_ae.s = self.bladegeom.s * sum(l) / self.bladegeom.s[-1]
        self.blade_ae.rthick = self.bladegeom.rthick * 100.
        self.blade_ae.chord = self.bladegeom.chord * self.blade_length
        self.blade_ae.aeset = np.ones(len(self.blade_ae.s))

    def calculate_c12axis(self):
        """
        compute the 1/2 chord axis based on the blade axis and chordwise
        rotation point

        nb: this only works for straight blades! #FIXME:
        """

        # The HAWC2 blade axis is defined using the 1/2 chord points
        b = self.bladegeom
        c12axis = np.zeros((b.x.shape[0], 4))

        for i in range(b.x.shape[0]):
            xc12 = (.5 - b.p_le[i]) * b.chord[i] * np.cos(b.rot_z[i] * np.pi / 180.)
            yc12 = - (.5 - b.p_le[i]) * b.chord[i] * np.sin(b.rot_z[i] * np.pi / 180.)
            c12axis[i, 0] = -(b.x[i] + xc12)
            c12axis[i, 1] = b.y[i] + yc12
            c12axis[i, 2] = b.z[i]
        c12axis[:, 3] = b.rot_z
        return c12axis


class HAWC2BeamStructureIDO(object):
    """
    Component to connect parameters when one of the parameter in the VarTree
    is an optimization variable (cannot connect twice the same parameter)
    """

    def execute(self):
        bps = HAWC2BeamStructure()
        self.h2beam_structure = []
        for name in bps.list_vars():
            fused_name = name
            if name == 'K': fused_name = 'J'
            try:
                val = getattr(self.beam_structure, fused_name)
                setattr(bps, name, val)
            except:
                pass
        self.h2beam_structure.append(bps)

if __name__ == '__main__':

    pass
