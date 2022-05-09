# Changelog

## 2022-05-04
Version 1.7.1
- added possibility to extrapolate aponeuroses from 100% or 50% of detected length.
- added fascicle length calculated from mean thickness and angle ("Fascicle_length_trig").

NB: The aponeurosis angle is the slope of the straight line fitted to the chosen portion of the aponeurosis (i.e. 50 or 100%). **It will influence the obtained pennation angle and fascicle length with any method**.

## 2020-03-17
Version 1.7
- added error handling.
- check depencies.
- fixed bug that caused analysis to fail when images are flipped and cropping manual.
