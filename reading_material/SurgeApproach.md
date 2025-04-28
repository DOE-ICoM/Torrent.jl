
## Modeling Storm Surge

### Determining Source Locations


### Driving Surge Inland

The location of a point $(x,y)$ in a new coordinate system rotated an angle $\theta$ about a point, $(x_0, y_0)$, in the original frame can be determined from the following formulas:
$$
  x^\prime = (x-x_0) \, cos(\theta) + (y-y_0) \, sin(\theta) \\
  y^\prime = -(x-x_0) \, sin(\theta) + (y-y_0) \, cos(\theta)
$$


### Surface Roughness

NLCD to surface roughness mapping available on page 4-76 of Hazus Hurricane Model Technical Manual, v5.1

could use the local roughness value, $z_0$, or could use an effective roughness value, $z^{eff}_0$, that results from some kind of averaging of $z_0$ along a, potentially 2D, fetch upwind of the actual grid cell.