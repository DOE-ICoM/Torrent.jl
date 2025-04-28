# Dam Failure Reservoir Dynamics

## Rectangular Cross Section

The following is based on the equation found [here](https://learn.hydrologystudio.com/studio-express/knowledge-base/weirs/).

Let's start with the basic Weir Equation for a rectangular cross section, where $Q$ is the flux over the weir and $h$ is the height of the water surface above the top of the weir.

$$
Q = C_w L h^{3/2}
$$

The height, $h$, is related to the total, time-dependent reservoir volume, $V$, remaining behind the weir. Let's assume that the reservoir has a rectangular bathymetry with a water depth, $z$, and surface area, $A$; in this case the reservoir volume is given by $V = A z$. Let's let the weir or dam height itself be given by, $z_0$, such that $h = z - z_0$. Since $z = V/A$, we have $h = V/A - z_0$.

And, of course, the outflow over the dam is just the change in volume of the water in the reservoir: $Q = -dV/dt$. Substituting these:
$$
\frac{dV}{dt} = - C_w L (V/A - z_0)^{3/2}
$$

## Trapazoidal Cross Section

Flux is based on the equation for $Q$ above, plus twice the following:
$$
Q^\prime = \frac{2}{5} C_w Z h^{5/4} \\
\frac{dV^\prime}{dt} = - C_w Z (V/A - z_0)^{5/4}
$$
where $Z$ is the side slope (horizontal to vertical) ratio.
