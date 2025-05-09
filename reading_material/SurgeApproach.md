
## Modeling Storm Surge

### Determining Source Locations

$$
\frac{dV}{dt} = \phi^{(top)} - \phi^{(side)} \\
l^2 \frac{dh}{dt} = \phi^{(top)} - h l U
$$

where $V$ is the volume of water in the cell, $\phi^{(top)}$ is the incoming flux introduced at the top of the cell, $\phi^{(side)}$, is the flux of water leaving out the side of the cell, $l$ is the length of a cell side, $h$ is the water depth, and $U$ is the fluid velocity out the cell side. Let's approximate $U$ using the Manning formula:
$$
U = \frac{1}{n} h^{2/3} \left( \frac{dh}{dx} \right)^{1/2}
$$
We then have:
$$
l^2 \frac{dh}{dt} = \phi^{(top)} - \frac{l h^{5/3}}{n} \left( \frac{dh}{dx} \right)^{1/2}
$$
Or flipping this around:
$$
\phi^{(top)} = l^2 \frac{dh}{dt} + \frac{l h^{5/3}}{n} \left( \frac{dh}{dx} \right)^{1/2}
$$
Perhaps we can find an anlytical expression for $dh/dx$ then:
$$
\frac{dh}{dx} = \frac{dh}{dt} \frac{dt}{dx} = \frac{dh}{dt} U^{-1} \\
\frac{dh}{dx} = \frac{dh}{dt} \frac{n}{h^{2/3}} \left( \frac{dh}{dx} \right)^{-1/2} \\
\frac{dh}{dx} = \left( \frac{dh}{dt} \frac{n}{h^{2/3}} \right)^{2/3}
$$
This then can be substituted back in, providing us with:
$$
\phi^{(top)} = l^2 \frac{dh}{dt} + \frac{l h^{5/3}}{n} \left( \frac{dh}{dt} \frac{n}{h^{2/3}} \right)^{1/3} \\
\phi^{(top)} = l^2 \frac{dh}{dt} + \frac{l \, h^{13/9}}{n^{2/3}} \left( \frac{dh}{dt}  \right)^{1/3}
$$
(Believe it or not, I think the units of that all work out.)

<!-- \phi^{(top)} = \left( l^2 + \frac{lh^{7/3}}{n^2} \right) \frac{dh}{dt} -->

But if we look at some storm surge depth observations [REF???] they tend of have a fairly Gaussian shape, that is the depth we want to affect has a functional form that looks something like:
$$
h(t) = h^* \exp \left( -\frac{(t-t^\prime)^2}{2 \, \sigma^2} \right)
$$
And, taking the derivative of this, we find:
$$
\frac{dh}{dt} = \frac{h^* (t^\prime-t) }{\sigma ^2} \exp \left(-\frac{(t-t^\prime)^2}{2 \sigma ^2} \right) \\[8pt]
\frac{dh}{dt} = \left( \frac{t^\prime-t}{\sigma^2} \right) h
$$
So, we have:
$$
\phi^{(top)} = l^2 \left( \frac{t^\prime-t}{\sigma^2} \right) h + \frac{l \, h^{16/9}}{n^{2/3}} \left( \frac{t^\prime-t}{\sigma^2} \right)^{1/3}
$$
Substituting the exponential form for $h$ into our equation for $\phi^{(top)}$, we have:
$$
\phi^{(top)} = l^2 \left( \frac{t^\prime-t}{\sigma^2} \right) h^* \exp \left( -\frac{(t-t^\prime)^2}{2 \, \sigma^2} \right)+ \frac{l}{n^{2/3}} \left( \frac{t^\prime-t}{\sigma^2} \right)^{1/3} \left[ h^* \exp \left( -\frac{(t-t^\prime)^2}{2 \, \sigma^2} \right) \right]^{16/9}
$$

<!--
Substituting these two back into the equation for $\phi^{(top)}$, we have:
$$
\phi^{(top)} = \left( l^2 + \frac{l h^{7/3}}{n^2} \right) \frac{h^* (t^\prime-t) }{\sigma ^2} \exp \left(-\frac{(t-t^\prime)^2}{2 \sigma ^2} \right) \\

\phi^{(top)} = \left( l^2 + \frac{l}{n^2} \left[ h^* \exp \left( -\frac{(t-t^\prime)^2}{2 \, \sigma^2} \right) \right]^{7/3} \right) \frac{h^* (t^\prime-t) }{\sigma ^2} \exp \left(-\frac{(t-t^\prime)^2}{2 \sigma ^2} \right)
$$
-->


### Driving Surge Inland

The location of a point $(x,y)$ in a new coordinate system rotated an angle $\theta$ about a point, $(x_0, y_0)$, in the original frame can be determined from the following formulas:
$$
  x^\prime = (x-x_0) \, cos(\theta) + (y-y_0) \, sin(\theta) \\
  y^\prime = -(x-x_0) \, sin(\theta) + (y-y_0) \, cos(\theta)
$$


### Surface Roughness

NLCD to surface roughness mapping available on page 4-76 of Hazus Hurricane Model Technical Manual, v5.1

could use the local roughness value, $z_0$, or could use an effective roughness value, $z^{eff}_0$, that results from some kind of averaging of $z_0$ along a, potentially 2D, fetch upwind of the actual grid cell.