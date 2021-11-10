# emacs: -*- mode: python-mode; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the NiBabel package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
"""
Geometric transforms.

.. currentmodule:: nitransforms

.. autosummary::
   :toctree: ../generated

   transform
"""
from .linear import Affine
from .nonlinear import DisplacementsFieldTransform
from .io import afni, fsl, itk


__all__ = ['afni', 'fsl', 'itk', 'Affine', 'DisplacementsFieldTransform']
