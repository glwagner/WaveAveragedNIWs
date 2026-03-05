# Gravitational, Symmetric, and Baroclinic Instability of the Ocean Mixed Layer

**Thomas W. N. Haine and John Marshall**

*Journal of Physical Oceanography*, Vol. 28, April 1998, pp. 634–658

---

<\!-- Page 1 -->
634
VOLUME 28
J O U R N A L O F P H Y S I C A L O C E A N O G R A P H Y
q 1998 American Meteorological Society
Gravitational, Symmetric, and Baroclinic Instability of the Ocean Mixed Layer
THOMAS W. N. HAINE* AND JOHN MARSHALL
Department of Earth, Atmospheric, and Planetary Sciences, Massachusetts Institute of Technology, Cambridge, Massachusetts
(Manuscript received 29 April 1996, in ﬁnal form 23 July 1997)
ABSTRACT
A hierarchy of hydrodynamical instabilities controlling the transfer of buoyancy through the oceanic mixed layer
is reviewed. If a resting ocean of horizontally uniform stratiﬁcation is subject to spatially uniform buoyancy loss
at the sea surface, then gravitational instability ensues in which buoyancy is drawn from depth by upright convection.
But if spatial inhomogeneities in the ambient stratiﬁcation or the forcing are present (as always exist in nature),
then horizontal density gradients will be induced and, within a rotation period, horizontal currents in thermal-wind
balance with those gradients will be set up within the mixed layer. There are two important consequences on the
convective process:
1) Upright convection will become modiﬁed by the presence of the thermal wind shear; ﬂuid parcels are exchanged
not along vertical paths but, rather, along slanting paths in symmetric instability. Theoretical considerations
suggest that this slantwise convection sets the potential vorticity of the mixed layer ﬂuid to zero but, in general,
will leave it stably stratiﬁed in the vertical.
2) The convective process ultimately gives way to a baroclinic instability of the horizontal mixed layer density
gradients. The resulting baroclinic waves are important agents of buoyancy transport through the mixed layer
and can be so efﬁcient that the convective process all but ceases.
The authors illustrate and quantify these ideas by numerical experiment with a highly resolved nonhydrostatic
Navier–Stokes model. Uniform spatial cooling at the surface of a resting, stratiﬁed ﬂuid in a 2½-dimensional model
on an f plane, in which zonal strips of ﬂuid conserve their absolute momentum, causes energetic vertical overturning.
A well-mixed boundary layer develops over a depth that is accurately predicted by a simple 1D law. In contrast,
differential surface cooling induces a mixed layer front. Fluid parcels, made dense at the surface, sink along slanting
trajectories in intense nonhydrostatic plumes. After cooling ceases the Ertel potential vorticity within the convective
layer is indeed found to be vanishingly small, corresponding to convective neutrality measured in the absolute
momentum surfaces that are tilted from the vertical by the horizontal vorticity of the thermal wind.
In analogous fully three-dimensional calculations, the absolute momentum constraint is broken, and the convection
at ﬁrst coexists with, but is ultimately dominated by, a baroclinic instability of the mixed layer. For typical mixed
layer depths of 500 m stability analysis predicts, and our explicit calculations conﬁrm, that baroclinic waves with
length scales O(5 km) develop with timescales of a day or so. By diagnosis of fully developed mixed layer turbulence,
the authors assess the importance of the baroclinic eddy ﬁeld as an agency of lateral and vertical buoyancy ﬂux
through the layer. A novel scaling for the lateral buoyancy ﬂux due to the baroclinic eddies is suggested. These
ideas are based on analysis of several experiments in which the initial stratiﬁcation, rotation rate, and buoyancy
forcing are varied, and the results are compared to previous attempts to parameterize the effects of baroclinic
instability. There is a marked difference between the scaling that accounts for the resolved experiments and the
Fickian schemes used traditionally in large-scale ocean models.
Finally, consideration of the results in light of high-resolution mixed layer hydrographic surveys in the northeast
Atlantic suggests mixed layer baroclinic instability may be very important at fronts. The authors speculate that the
process exerts a large inﬂuence on the character of newly subducted thermocline water throughout the extratropical
ocean.
1. Introduction
Knowledge and accurate representation of the pro-
cesses controlling the development of the upper ocean
* Current afﬁliation: Atmospheric, Oceanic and Planetary Physics,
Clarendon Laboratory, Oxford, United Kingdom.
Corresponding author address: Dr. Thomas W. N. Haine, Atmo-
spheric, Oceanic and Planetary Physics, Clarendon Laboratory, Parks
Road, Oxford OX1 3PU, United Kingdom.
E-mail: T.Haine1@physics.ox.ac.uk; marshall@gulf.mit.edu
is vital if we are to understand the ways in which the
large-scale ocean structure is determined and develop
a quantitative theory of it. The surface mixed layer of
the ocean, directly in contact with the atmosphere, is
of central importance in determining the manner and
rates of heat, freshwater, momentum, and gas exchange
with the interior of the ocean. Traditional oceanic
mixed layer paradigms (e.g., Kraus and Turner 1967;
Mellor and Yamada 1974) suppose that the properties
of this surface layer are set by vertical mixing caused
by mechanical stirring from the wind, surface gravity
wave breaking and by convective mechanisms induced
Unauthenticated | Downloaded 03/05/26 01:04 AM UTC


<\!-- Page 2 -->
APRIL 1998
635
H A I N E A N D M A R S H A L L
FIG. 1. (a) A traditional mixed layer model where a layer with zero vertical stratiﬁcation is developed by
gravitational overturning and surface buoyancy loss. (b) In the presence of a lateral density gradient in a
rotating frame, symmetric instability sets the potential vorticity to zero, leaving weak stratiﬁcation. This layer
is unstable to baroclinic waves, which results in lateral buoyancy transfer in the mixed layer.
by buoyancy loss from the sea surface. The latter pro-
cess dominates in mixed layers deeper than O(100 m)
and is a consequence of the familiar gravitational over-
turning (or upright convection) associated with denser
ﬂuid overlying lighter ﬂuid (represented schematically
in Fig. 1a).
The upper ocean is not horizontally homogeneous,
however, as is clearly revealed by any high-resolution
survey (e.g., Samelson and Paulson 1988). For ex-
ample, Fig. 2 is a section from the northeast Atlantic
(obtained in April 1991 by a SeaSoar—a towed, un-
dulating CTD; Cunningham et al. 1992; Pollard 1986)
revealing density gradients from the smallest resolved
scales [O(5 km) horizontally: O(5 m) vertically] to the
length and depth of the survey. One frequently ob-
serves unstable regions adjacent to stratiﬁed ﬂuid. At
the northern end of the section the upper 250 m is
homogeneous, or slightly unstable, a signature that ac-
tive overturning is under way. A warmer layer is seen
to the south that has large regions of very weak strat-
iﬁcation despite a surface cap of lighter ﬂuid in the
upper 20–30 m. At 45.48N, near the southern end of
the section, an anticyclonic eddy is apparent. This fea-
ture occupies the upper 300 m of the water column
and has a characteristic diameter of 25 km. The core,
and the surrounding ﬂuid, is weakly stratiﬁed with stat-
ically unstable patches. A shallower (upper 100 m)
feature is present at 46.58N and is almost detached from
the less dense water to the south. These observations
give no obvious indication of a vertically homogeneous
mixed layer separated from stratiﬁed water below. In-
deed, it is very difﬁcult to deﬁne the mixed layer in
an unambiguous way. The mixed layer depth, diag-
nosed as the depth at which the density exceeds the
surface value by 0.05 kg m23, is shown, but does not
correspond to any clear mixed layer base. In fact, the
mixed layer depth determined in this way is very sen-
sitive to the exact criterion used. Nevertheless, there
are signiﬁcant lateral gradients within the convectively
stirred layer, caused by a variety of hydrodynamical
processes induced by surface buoyancy and momentum
ﬂuxes. Clearly, if these hydrodynamical processes are
sufﬁciently slow and large scale, the earth’s rotation
will inﬂuence them.
In this paper we review and investigate some of the
key processes that control the ﬂux of buoyancy ver-
tically and horizontally in the upper ocean. We argue,
and illustrate by numerical experiment, that in the pres-
ence of lateral density gradients upright convection can
be modiﬁed by thermal wind shear so that overturning
occurs along paths that slant to the vertical. This slant-
wise convection rapidly (typically over a few hours)
restores the Ertel potential vorticity of the convecting
layer to zero and maintains a layer with weak vertical
stratiﬁcation. But of equal importance is that this state
is susceptible to nonhydrostatic baroclinic instability,
which quickly develops causing vertical and lateral
transfer of buoyancy in the mixed layer (Fig. 1b) on
geostrophic scales. Numerical experiments show that
baroclinic instability in the mixed layer can result in
lateral buoyancy ﬂuxes that signiﬁcantly modify the
shoaling and deepening of the layer and are so efﬁcient
that the convective process all but vanishes.
Finally, scaling laws are deduced, and tested against
Unauthenticated | Downloaded 03/05/26 01:04 AM UTC


<\!-- Page 3 -->
636
VOLUME 28
J O U R N A L O F P H Y S I C A L O C E A N O G R A P H Y
FIG. 2. Vertical section of potential density (kg m23: referenced to the surface) along a southward cruise track in the northeast Atlantic
in April 1991 (Vivaldi cruise). The measurements were made from a SeaSoar instrument (Pollard 1986), a towed undulating CTD that cycles
between the surface and 500 m every 4 km along the track. The dashed line shows the mixed layer depth diagnosed as the depth at which
the density exceeds the surface density by 0.05 kg m23. The contours mark those regions where the vertical stratiﬁcation is very weak or
unstable [N 2 , 1 3 1028 s22, where N is the Brunt–Va¨isa¨la¨ frequency, Eq. (2)]. The proﬁles of N 2 were smoothed with a 20-db low-pass
ﬁlter.
our numerical experiments, in an attempt to quantify
the importance of these dynamical processes in the
ﬁeld. They suggest that lateral transfer by mixed layer
baroclinic instability may be ubiquitous and represent
an important process that is absent from one-dimen-
sional mixed layer paradigms. It is most signiﬁcant at
the time of subduction, in the aftermath of deep-reach-
ing convection (shortly after the data in Fig. 2 were
taken), and may be inﬂuential in setting the charac-
teristics of newly formed thermocline water.
2. Theoretical background
The various types of instability that occur naturally
in the oceanic mixed layer are hybrid in nature. It is
instructive, however, to consider idealized limit cases
in which one or another of the destabilizing forces acts
alone. Four simple types of instability are of interest
(they are reviewed in a meteorological context by Eady
1951) and are discussed here in the context of oceanic
mixed layers:
1a) Gravitational instability (‘‘ordinary’’ convection
or static instability)
1b) Centrifugal instability (sometimes called inertial
instability)
2a) Baroclinic instability (of the thermal wind)
2b) Helmholtz instability (at a velocity discontinuity).
Instability of type 1a has been the traditional focus in
studies of the mixed layers. Instability of type 1b is un-
likely to occur in pure form in the mixed layer, but the
hybrid gravitational–centrifugal instability 1ab (known as
‘‘symmetric instability’’) is likely to occur and is one of
the focuses of the present study. Instability of type 2a is
not widely associated with mixed layer processes, but a
number of theoretical studies (Stone 1971; Young and
Chen 1995; Fukamachi et al. 1995; Barth 1994) have
pointed to its importance. We argue here that baroclinic
instability may be one of the most important processes
active in mixed layers. Finally, instability of type 2b is
signiﬁcant in mixed layers where there is a strong ver-
tically sheared ﬂow—especially at low latitudes (see Pol-
lard et al. 1973; Price et al. 1986). Because it is a me-
chanically driven, rather than a buoyancy driven, process
we do not consider it further here. In the following brief
review we make much use of parcel theory—this is sum-
marized in appendixes A and B.
a. Type 1a: Gravitational instability
Consider a resting ocean of constant stratiﬁcation Nth
(the thermocline stratiﬁcation) subject to uniform and
widespread buoyancy loss from its upper surface as
shown in Fig. 1a. The thermodynamic equation is
Unauthenticated | Downloaded 03/05/26 01:04 AM UTC


<\!-- Page 4 -->
APRIL 1998
637
H A I N E A N D M A R S H A L L
Db 5 B,
(1)
Dt
where b 5 2(g/ro)s is the buoyancy; s the potential
density; r0 a constant reference density; g the acceler-
ation due to gravity; B 5 ]B/]z is the buoyancy forcing,
the divergence of the ﬂux B; and D/Dt is the material
derivative.
The ﬂuid cannot simultaneously overturn on the large
scale; rather the qualitative description must be that the
response to widespread cooling is one in which rela-
tively small convection cells (plumes) develop. Fluid
parcels at the surface become dense and sink under grav-
ity displacing less dense parcels from below. The con-
tinual exchange of ﬂuid parcels in this way will, over
time, create a layer that, as we show below, is very
close to neutral with respect to its thermodynamic prop-
erties. However, as long as the buoyancy loss persists
there will be, on the average, a small statically unstable
buoyancy gradient:
]b
2
N
5
, 0,
(2)
mix
]z
where Nmix is the Brunt–Va¨isa¨la¨ frequency in the over-
turning layer. If
, 0, then exchange of parcels
2
Nmix
vertically must release gravitational potential energy.
Since horizontal motion does not affect potential energy,
one need only consider vertical overturning; parcel the-
ory [see appendix A, Eq. (A3)] then yields an upper
limit on the growth rate v:
v2 # |
|.
2
Nmix
(3)
For a prescribed
, with appropriate boundary con-
2
Nmix
ditions, it is straightforward to obtain complete solutions
through linear stability analysis (see, e.g., Rayleigh
1916; Veronis 1958; Chandrasekhar 1961); these show
that v2 is nearly attained when the convection cells are
tall and thin, for which little energy is supplied to hor-
izontal motion. Laboratory simulations, however, sug-
gest that the aspect ratio of fully developed turbulent
convection approaches unity such that horizontal and
vertical scales are of the same order.1
Many competing effects combine to control the de-
tailed dynamics on the plume scale. However, irrespec-
tive of these details, the gross transfer properties of the
population of convective cells must be controlled by the
large scale; the raison d’eˆtre for the overturning is to
ﬂux buoyancy vertically to offset buoyancy loss at the
surface. As shown in appendix A, the following ‘‘law’’
1 The detailed dynamics setting the plume scale in ocean convection
is, as yet, unclear. For example, if the convecting layer becomes deep
enough then the concomitant increase in the lateral scale may cause
the convection to be inﬂuenced by Taylor columns associated with
the earth’s rotation. This may occur in open-ocean deep convection,
but will not be pursued here. In recent years it has been studied at
some length—see the review by Marshall et. al. (1994).
of vertical buoyancy transfer for the plume scale Bp can
be developed using the same parcel theory that leads to
Eq. (3):
Bp 5 wDb 5 Dz1/2Db3/2,
(4)
where w is the vertical velocity in the plume, Dz its
vertical extent, and Db is the difference in buoyancy of
the rising and sinking ﬂuid. Note that the above scaling
is appropriate in the highly supercritical limit and as-
sumes that visco–diffusive parameters are irrelevant. It
also assumes that the mixed layer is sufﬁciently shallow
that rotational effects are not dominant; that is, Ro* 5
(1/h)
B/ f 3, the natural Rossby number, is sufﬁciently
Ï
large (h is the mixed layer depth and f the Coriolis
parameter). Then the classical nonrotating scalings are
appropriate. In some of the convection experiments of
Jones and Marshall (1993) Ro* was small enough for
rotation to affect the convection scale. In the experi-
ments presented below the rotational limit does not ap-
ply.
If the plumes, acting in concert, achieve a vertical
buoyancy ﬂux sufﬁcient to balance loss from the sur-
face, then Bp 5 B0, and choosing a deep mixed layer
exposed to a heat loss of ;500 W m22, typical of deep
winter mixed layers:
Dz 5 1000 m;
B0 5 1027 m2 s23,
we deduce that the temperature anomalies, DT ;
0.0018C, typical vertical velocities are [using (A2) of
the appendix A] a few centimeters per second with
timescales of perhaps 8 h or so. Here we have assumed
that the buoyancy loss is all due to heat and the thermal
expansion coefﬁcient of water a is 2 3 1024 K21. It is
notable that such a tiny temperature difference between
rising and sinking ﬂuid parcels can achieve a very large
heat (and buoyancy) ﬂux. We conclude that in the ab-
sence of lateral inhomogeneities the vertical column
within the convecting layer is indeed very well mixed;
in Eq. (2) is very small relative to typical ther-
2
Nmix
mocline stratiﬁcations. In the limit that Nmix/Nth K 1
and to the extent that entrainment of stratiﬁed ﬂuid from
the base of the mixed layer can be neglected, Eq. (1)
tells us that the depth of the mixed layer h must increase
with time t according to (Turner 1973):
Ï2B t
0
h 5
,
(5)
Nth
where B0 is the buoyancy forcing at the sea surface and
Nth is the stratiﬁcation of the underlying ﬂuid.
b. Type 1b: Centrifugal instability
Suppose that the mixed layer is of uniform density
everywhere (
5 0) and that a horizontal barotropic
2
Nmix
velocity u(y) exists within the layer with du/dy a con-
stant. Now energy exists only in kinetic form. Let us
consider overturning in the yz plane with no variations
Unauthenticated | Downloaded 03/05/26 01:04 AM UTC


<\!-- Page 5 -->
638
VOLUME 28
J O U R N A L O F P H Y S I C A L O C E A N O G R A P H Y
in x. Then a consideration of the zonal momentum equa-
tion (assumed inviscid) on an f plane tells us that
Dm 5 0,
(6)
Dt
where
m 5 u 2 f y
(7)
is the absolute momentum with f the Coriolis parameter
and u the zonal velocity.
Since u 5 u(y), the m surfaces are vertical at all times
and parcel theory (appendix B) tells us that energy is
released if uy . f, corresponding to negative absolute
vertical vorticity z, and that the growth rate is
v2 # f(uy 2 f) 5 2 fz.
(8)
Linear theory predicts that the maximum growth rate
occurs for shallow, broad cells since there is no energetic
advantage to be gained in exchanging particles verti-
cally.
Pure centrifugal instability of this kind is not likely
to be common in oceanic mixed layers because of the
stabilizing effect of the earth’s rotation, except perhaps
at low latitudes in regions of anticyclonic vorticity.
c. Type 1ab: Symmetric instability;
gravitational–centrifugal instability
Now let us combine the results of sections 1a and 1b
to consider a system in which neither
or uy vanishes.
2
Nmix
This corresponds to the normal state of affairs in a mixed
layer—drawn schematically in Fig. 1b—in which the
density varies in the horizontal across the mixed layer
(because of more vigorous convection on one side than
the other, for example). On the large-scale a zonal cur-
rent u(y, z) will develop in thermal wind balance with
lateral mixed layer density gradients, given by
by
u 5 m 5 2
.
(9)
z
z
f
The presence of rotation and a zonal ﬂow in thermal
wind balance with a lateral density gradient place,
through Eq. (6), important rotational and angular mo-
mentum constraints on the convective process. The m
surfaces, which before were vertical, are now, because
of the presence of the thermal wind and its associated
horizontal component of vorticity, tilted over. Because
the m surfaces are material surfaces, they will induce
ﬂuid particles to move along slanting rather then vertical
paths. Moreover, the stability of the layer will depend
on the sign of =b measured in the m surface (corre-
sponding to gravitational instability) or the sign of the
absolute vorticity normal to the b surface (correspond-
ing to centrifugal instability). Both viewpoints are com-
plementary and entirely equivalent. Emanuel (1994)
calls this more general mixed instability ‘‘slantwise con-
vection.’’ The stability depends on the sign of the po-
tential vorticity (Hoskins 1974):
1
Q 5
h·=b,
(10)
g
a measure of the stratiﬁcation in the direction of h, the
absolute vorticity vector or, equivalently, a measure of
h normal to b surfaces. If, as in our thought experiment,
there are no variations in x, then the absolute vorticity
vector lies in a surface of constant absolute momentum,
and Q is just the Jacobian of m and b:
1
Q 5
J (m, b).
(11)
yz
g
If Q is negative, then the ﬂow is unstable to symmetric
instabilities and slantwise convection might be expected
to return the Q of the layer to zero, the state of marginal
stability. The sign of Q depends on the slope of the m
surfaces relative to the b surfaces and is zero when they
are exactly coincident; in the limit of zero Q there is
no stratiﬁcation in an m surface and the component of
h normal to the b surface is zero. The magnitude of the
absolute vorticity, resolved perpendicular to the b sur-
faces, is simply |h| 5 gQ/|=b|. For small slopes |=b| ;
and
2
Nmix
2
f N
z
1
mix
Q 5
2
,
(12)
1
2
g
f
Ri
where z 5 f 2 uy is the vertical component of the
absolute vorticity and Ri 5
/
is the Richardson
2
2
N
u
mix
z
number.
Parcel theory can be readily employed to analyze the
stability of a zonal ﬂow in thermal wind balance to
overturning in a vertical plane. The method is outlined
in appendix B—see also chapter 12 of Emanuel (1994).
Maximum release of energy occurs when ﬂuid parcels
are exchanged along surfaces coincident with the b sur-
faces. Then parcel theory yields
z
1
2
2
v # 2 f
2
,
(13)
1
2
f
Ri
and the ﬂow will be unstable when Q , 0 or, equiva-
lently,
Ri , f/z
(14)
[see Eq. (B5), appendix B].
Thus, we see that symmetric instability (overturning
in the vertical plane with conservation of zonal mo-
mentum) is to be expected when
R absolute vorticity is small (anticyclonic shear), though
it need not be negative
R horizontal thermal gradient is strong (large uz)
R the static stability is small.
These conditions are likely to be met frequently within
the oceanic mixed layer.
Moreover, the arguments above suggest that, if Q ,
0, then one might expect convection (appropriately gen-
eralized in the ‘‘symmetric’’ sense) to occur and that
Unauthenticated | Downloaded 03/05/26 01:04 AM UTC


<\!-- Page 6 -->
APRIL 1998
639
H A I N E A N D M A R S H A L L
FIG. 3. Dispersion relation for baroclinic instability in a uniformly
stratiﬁed layer where Ri 5 1. The dashed line is Eady’s (1949) result,
appropriate for asymptotic large Ri; the full line is due to Stone’s
(1971) nonhydrostatic theory for Ri $ 1, where the ratio of the aspect
ratio to the Rossby number is unity. The wavenumber has been non-
dimensionalized with a characteristic scale u/ f where u is the char-
acteristic zonal speed. The growth rate has been nondimensionalized
by the inertial frequency and the dimensional dispersion curve is for
the case of a layer 500 m deep.
the end state of the convective process will be one in
which Q →0.
d. Baroclinic instability
Let us now allow variations in x; the momentum of
ﬂuid parcels can be changed by pressure gradient forces
and will not be conserved. On the large-scale Coriolis
effects are important and, through the thermal wind re-
lation, make possible a storage of potential energy on
which instabilities can feed. The three spatial dimen-
sions permit currents in the x direction that are side-by-
side in a baroclinic instability, rather than above one
another as in the case of uniform overturning in sym-
metric instability. If the Richardson number of the ﬂow
is large, then the ensuing motion is quasigeostrophic;
the most unstable modes, on the scale of the Rossby
radius of deformation, exchange parcels along surfaces
that have a slope one-half of that of the isentropic sur-
face. The growth rate can be deduced using parcel the-
ory [see (B8)]2:
2f
2
v #
,
if
Ri k 1.
(15)
Ri
In the oceanic mixed layer Ri will not be large—it
is likely to be of order unity—and the quasigeostrophic
result (15) must be modiﬁed.
Stone (1971, 1972) derived expressions for the
growth of linear baroclinic waves in the low Richardson
number limit and showed that the results of Eady (1949)
apply with only quantitative modiﬁcation [essentially
Ri →(1 1 Ri)]. The baroclinic instability mechanism
endures in ﬂows in which Ri is O(1) and, indeed, co-
exists with symmetric instability if Ri , 1. Stone
showed that the growth rate of the most unstable (mod-
iﬁed) Eady mode is
f
v 5 0.304
,
(16)
1/2
(1 1 Ri)
and the scale of this mode is
Nh
1/2
L 5 1.016
(1 1 Ri)
Ri $ 1.
r
f
Figure 3 compares the Eady (1949) and Stone (1971)
dispersion relations for unstable modes when Ri 5 1.
2 Care must be taken in the application of parcel theory to baroclinic
instability because ﬂuid parcels do not conserve momentum, so com-
putations of kinetic energy change cannot be made with precision.
Moreover, unlike convection and symmetric instability, which are
local in nature, baroclinic instability is global and intimately asso-
ciated with the temperature distribution along the boundary. Nev-
ertheless analogies with convection can be usefully drawn. Indeed,
Eady (1949) drew such analogies and, along with exact solutions of
the linear problem, also showed how to deduce Eq. (15) heuristically
using parcel theory.
Although Eady’s analysis is formally inapplicable, his
result is qualitatively correct. Stone’s more general the-
ory predicts that the growth rates are slower and the
fastest growing modes are at larger scales. Figure 3 also
shows the corresponding dimensional values for a layer
of depth 500 m. In such an unstable mixed layer the
linear theory predicts that the fastest growing mode has
a wavelength of a few kilometers and will grow ex-
ponentially with a characteristic timescale of around a
day.
Summary. We conclude that there is a hierarchy of
instability mechanisms potentially at work in oceanic
mixed layers; theory suggests that symmetric instability
and baroclinic instability ought to be ubiquitous in
mixed layers and potentially important in the mixed
layer buoyancy budget. In the following sections we
present numerical calculations that conﬁrm the impor-
tance of these processes and illustrate
1) the evolution of a symmetric instability when ab-
solute momentum is conserved by zonal strips of
ﬂuid, showing convection along slanting paths and
the restoration of the mixed layer to the neutral, zero
potential vorticity state
2) the development of a baroclinic instability of the zero
potential vorticity state when the absolute momen-
tum constraint is relaxed
3) the modiﬁcation of the properties of the mixed layer
due to lateral and vertical eddy transfer of buoyancy
due to fully developed geostrophic eddies.
Unauthenticated | Downloaded 03/05/26 01:04 AM UTC


<\!-- Page 7 -->
640
VOLUME 28
J O U R N A L O F P H Y S I C A L O C E A N O G R A P H Y
3. Numerical experiments in convective,
symmetric, and baroclinic instability
In order to test the theoretical ideas summarized in
the previous section we employ a numerical model con-
ﬁgured to focus on the processes of interest. The model,
described in Marshall et al. (1997a,b) solves the Navier–
Stokes equations for a Boussinesq incompressible ﬂuid
using ﬁnite-volume techniques; it need not make the
hydrostatic approximation.
The model geometry is shown in Fig. 4a. A periodic
channel with constant depth is used of, nominally, length
50 km and width 30 km. The cell dimension is 250 m
in the horizontal with vertical spacing varying between
40 m at the surface and 400 m at the bottom. Initially
the ﬂuid is uniformly stratiﬁed and motionless. We adopt
a linearized equation of state with one thermodynami-
cally active variable:
s 5 r0[1 2 a(T 2 T0)],
(17)
where the expansion coefﬁcent a is 2 3 1024 K21 at
temperature T0.
To represent unresolved dynamics and ensure nu-
merical stability a Laplacian diffusion of heat and
momentum is applied. The diffusivities and viscosi-
ties are equal with horizontal and vertical magnitudes
of 5 and 0.02 m 2 s21, respectively. The horizontal
diffusivity was chosen to be the smallest possible to
guarantee a coherent vertical velocity and vorticity
on the grid. Our numerical solutions are less sensitive
to the vertical diffusivity; indeed the vertical diffu-
sivity could have been set to zero without any dele-
terious effects. Free slip is allowed at the solid bound-
aries. The model also integrates a dynamically passive
tracer held at a constant value in the uppermost layer,
and initially set to zero elsewhere. The calculations
presented in this paper were carried out on a 128-
node CM5 computer at the Massachusetts Institute of
Technology.
a. Gravitational instability; upright convection
First, we investigate upright overturning using 2D (y,
z) nonhydrostatic dynamics where initially the ﬂuid is
uniformly stratiﬁed (Nth 5 8 3 1024 s21), resting, and
rotating at a constant rate of 1024 s21. The motion is
forced by a steady, constant, surface buoyancy loss of
2 3 1027 m2 s23, that corresponds to an oceanic heat
loss of 400 W m22 (Table 1: experiment 1).
Results after 9 days show energetic vertical over-
turning in a boundary layer several hundred meters thick
(Fig. 5a; upper 1000 m only), whereas below there is
extremely weak ﬂow—essentially this region is undis-
turbed. The gross aspects of the convection are resolved
in the model, albeit coarsely. The streamfunction shows
the convection cells are organized into rolls with no
preferred direction of slant in the vertical plane. Figure
5b is a horizontal average of the boundary layer in Fig.
5a. It reveals a temperature inversion close to the surface
of around 0.0258C, whereas in the interior of the con-
vective layer there is a much smaller vertical temper-
ature gradient with a contrast of 0.0038C over the full
depth of the boundary layer, in good accord with the
parcel theory of section 1a. The gross predictions from
the ‘‘law of vertical buoyancy ﬂux’’ [Eq. (4)] derived
in appendix A are well supported by the explicit cal-
culation.
Using the mean temperature proﬁle (Fig. 5b) we es-
timate a depth for the mixed layer using a criterion based
on the stratiﬁcation. This procedure is applied at several
times during the evolution of the experiment and the
resulting time series plotted in Fig. 5c along with the
prediction of the simple 1D, nonpenetrative law [the
solid curve—Eq. (5)]. We see good agreement with the
theory. The conﬁrmation of the 1D result suggests that,
at least in this numerical model, the mixed layer deepens
and cools without signiﬁcant entrainment of the under-
lying water or formation of a front at the base of the
layer (in other words, there is no penetrative convec-
tion).
b. Symmetric instability
Now we examine convection in the presence of lateral
inhomogeneities that induce lateral density gradients and
a thermal wind in balance with it. Moreover, we conﬁgure
the model so that zonal strips of ﬂuid conserve their
absolute momentum (all ]/]x terms are set to zero). A
buoyancy ﬂux at the sea surface is speciﬁed that varies
across the channel according to a hyperbolic tangent,
thus,
B 5 B1/2{tanh[2(y 2 Ly/2)/Lf] 1 1},
(18)
where B1/2 is the buoyancy ﬂux at midchannel, y the
distance across the channel, Ly the channel width, and Lf
a characteristic length scale of the forcing (see Fig. 4b
and Table 1: experiment 2). The tanh function smoothly
changes the forcing across the channel and provides a
well-deﬁned maximum gradient in ﬂux, located at the
channel center. This allows a mixed layer to grow that
is deeper on one side of the channel than on the other,
inducing a lateral density gradient and a thermal wind in
balance with it.
Figure 6 shows the ﬁelds from the central portion of
the channel after 9 days of cooling. It is clear from the
isotherms that the overturning motions cause ﬂuid to
move systematically in slanting paths and maintain a non-
vanishing stratiﬁcation in the region that is being actively
mixed. However, the temperature ﬁeld alone is a poor
indicator of the regions of active overturning. Rather,
potential vorticity is the key dynamical variable as shown
by Fig. 6b. There are distinct plumes, of negative PV,
draining the surface source of negative PV into the in-
terior. The tracer in Fig. 6c shows striking similarities to
the PV distribution reminding us that both quantities are
materially conserved.
Unauthenticated | Downloaded 03/05/26 01:04 AM UTC


<\!-- Page 8 -->
APRIL 1998
641
H A I N E A N D M A R S H A L L
FIG. 4. (a) Schematic diagram of the model domain with an iso-
pycnal outcropping into the mixed layer with a lateral density gra-
dient. The Boussinesq, incompressible Navier–Stokes equations are
solved starting from a resting ﬂuid with uniform stratiﬁcation. (b)
Hyperbolic tangent cooling function used to generate a mixed layer
front.
The timescale of the symmetric instability has been
estimated using parcel theory. The relevant expression
for the growth rate is Eq. (13), which we now write as
2
2
f Nth
2
v # 2
Q*,
(19a)
2
Nmix
where Q* is the PV nondimensionalized by its thermo-
cline value
Q* 5 gQ/(
)
2
fN th
(19b)
and is the quantity plotted in Fig. 6b. Typical values of
/
and Q* are 4 and 20.3 respectively, yielding a
2
2
N
N
th
mix
growth rate of around 3 h.
The integration shown in Fig. 6 was continued for
another 24 hours, but now with the surface cooling
switched off. Figure 7 shows that within this period
almost all of the convective, turbulent ﬂow has
ceased, leaving a layer with nonvanishing vertical
stratiﬁcation but very small potential vorticity
(around 1% of the undisturbed value). This is con-
ﬁrmed by the close alignment of the isotherms with
the contours of absolute momentum m. The plumes
of negative PV have been mixed away by the sym-
metric instability, erasing density gradients along ab-
solute momentum surfaces and setting the Richardson
number to unity.
Also shown in these ﬁgures is the predicted mixed
layer depth, h, based on the simple one-dimensional
scaling Eq. (5). Inspection of Figs. 6 and 7 clearly
shows that the prediction is in good agreement with
the base of the mixed layer. Although the mixed layer
ﬂuid has nonzero N, that N is much smaller than Nth,
Unauthenticated | Downloaded 03/05/26 01:04 AM UTC


<\!-- Page 9 -->
642
VOLUME 28
J O U R N A L O F P H Y S I C A L O C E A N O G R A P H Y
TABLE 1. The numerical experiments. Experiments 1 and 2 were in a 2D domain (no x variation), experiments 3–16 fully 3D, experiments
1–13 were in a channel 30 km wide: experiments 14–16 in a channel 60 km wide. The Burger number quoted [Eq. (27)] is at time t 5
tmodel.
Experi-
ment
Nth
(31024 s21)
f
(31024 s21)
Cooling
width, Lf
(km)
Heat ﬂux
(W m22)
Buoyancy
ﬂux, B1/2
(31027 m2 s23)
By
(310211 m s23)
Bu
(31023)
lrot
(m)
tmodel
(days)
1
2
3
4
5
6
8.37
8.37
8.37
8.37
8.37
16.7
1.0
1.0
1.0
1.0
1.0
1.0
`
10
10
10
10
10
400
400
400
800
200
400
1.96
1.96
1.96
3.92
0.981
1.96
0
3.92
3.92
7.84
1.96
3.92
—
—
82.4
153
48.0
73.7
443
443
443
626
313
443
—
—
9.72
9.04
11.3
8.69
7
8
9
10
11
8.37
4.18
8.37
8.37
8.37
1.0
1.0
1.0
2.0
0.5
10
10
10
10
10
100
400
50
200
200
0.491
1.96
0.245
0.981
0.981
0.981
3.92
0.491
1.96
1.96
37.7
82.4
28.6
17.9
160
221
443
157
111
886
17.8
9.72
27.0
16.9
9.43
12
13
14
15
16
16.7
8.37
8.37
8.37
8.37
1.0
2.0
1.0
1.0
1.0
10
10
20
20
10
50
400
400
200
400
0.245
1.96
1.96
0.981
1.96
0.491
3.92
1.96
0.981
3.92
33.8
24.9
35.1
24.2
84.3
157
157
443
313
443
31.9
11.8
16.6
22.8
9.95
so there is little systematic departure from the 1D
scaling—the slantwise convection has resulted in lit-
tle cross-channel buoyancy ﬂux—a matter we ex-
amine in detail in section 4.
c. Baroclinic instability
We now perform an identical experiment, but in a
3D domain, allowing zonal as well as meridional vari-
ations, and hence the possibility of baroclinic insta-
bility (Table 1: experiment 3). A typical example of
the ﬂow development is shown in Fig. 8. The near-
surface ﬁelds of temperature reveal progression from
plume-scale convection at day 3 to ﬁnite amplitude
baroclinic instability at day 6 with a mature ﬁeld of
geostrophic turbulence by day 9. A surface-intensiﬁed
jet evolves in balance with the across-channel tem-
perature gradient, with the eddying part of the ﬂow
dominating. Since there is no stress applied at the
ocean surface, the global zonal momentum cannot
change and eastward ﬂow at the surface is compen-
sated by a westward current below. The length scale
for the baroclinic instability at day 6 is around 5 km—
somewhat larger than the prediction of Stone’s linear
instability analysis of about 3 km (Fig. 3). We believe
this is primarily due to a nonlinear inverse cascade
to larger scales familiar in two-dimensional turbu-
lence (Rhines 1975; Held and Larichev 1996, and ref-
erences therein).3
3 It is possible that the dissipation in the model may also inﬂuence
our results, although the broad conclusions are independent of the
size of the assumed diffusivities provided that they are sufﬁciently
small. Stability theory suggests (e.g., Lin and Pierrehumbert 1988)
that the wavelength of the fastest growing mode is insensitive to the
assumed viscosity, although the growth rate does show some sensi-
Figure 9 shows the zonal-mean sections of temper-
ature and absolute momentum. At day 3 the one-di-
mensional prediction agrees reasonably well with the
mixed layer depth one might diagnose assuming a ver-
tically homogeneous mixed layer. Indeed, the observed
stratiﬁcation is very weak except near the surface where
a static instability prevails, triggering the convective
instability. By day 6, however, weak stratiﬁcation has
returned to the mixed layer despite the ongoing buoy-
ancy loss. This result cannot be explained using one-
dimensional ideas, and by day 9 it is unclear how to
distinguish between mixed layer ﬂuid and the under-
lying water based on the temperature ﬁeld alone.
In Fig. 10a we plot the PV at a depth of 671 m on
day 9 and can identify the major dynamical processes
at work. The largest-scale features are the baroclinic
eddies, which interleave high PV ambient ﬂuid, typ-
ical of the southern part of the domain, with the con-
vectively modiﬁed low PV water to the north. There
is a strong gradient between these two water types
with relatively small volumes of water with inter-
mediate PV (Fig. 10b shows a volumetric census of
the various PV and N 2 classes in Fig. 10a). To the
north the patches of negative PV identify the energetic
plumes, which draw the PV and buoyancy from the
ﬂuid. Figure 11, a north–south vertical section at this
time, shows the motion is dominated by geostrophic
eddies over most of the channel: the ﬂow is quasi-
two-dimensional and coherent over several kilome-
ters. In contrast we observe regions of active con-
vection toward the northern wall. These nonhydro-
tivity, being reduced at higher values. Nonlinear interactions rapidly
result in structures with larger scales than linear theory predicts,
which are thus less susceptible to damping.
Unauthenticated | Downloaded 03/05/26 01:04 AM UTC


<\!-- Page 10 -->
APRIL 1998
643
H A I N E A N D M A R S H A L L
FIG. 5. Results from the gravitational overturning experiment (experiment 1—Table 1). (a) Vertical section
at day 9 showing isotherms (solid, contour spacing of 0.028C), overturning streamfunction (dashed, contour
spacing of 2 m2 s21) and ﬂow. The ﬂow is shown by the small dashes, which indicate equivalent displacements
after 30 min. The peak speeds are (0.069, 0.024) m s21 in the (y, z) directions. The thick dashed line is the
prediction of the 1D law for the depth of the mixed layer [Eq. (5)]. Only the central third and upper 1000 m
of the 30-km-wide and 2000-m-deep channel are shown. (b) Across-channel mean temperature proﬁle at day
9. (c) Time series of mixed layer depth. Full line is 1D prediction [Eq. (5)]—circles are model results (a vertical
temperature gradient of 5 3 1025 8C m21 is used to diagnose the mixed layer depth).
static plumes are distinctly different with O(1) aspect
ratio, large vertical velocities, and negative PV. It is
hard to distinguish between symmetric instability and
upright overturning, unlike in the calculations of sec-
tion 3b which clearly show the slantwise motions. The
effect on PV is clear however: It is rapidly reset to
zero by the convection induced by surface buoyancy
loss.
One day after the cessation of cooling the plumes of
negative PV have disappeared, leaving a large pool of
ﬂuid with very small potential vorticity (Fig. 12). This
low PV ﬂuid is then vigorously folded in with the un-
modiﬁed water to the south by eddies that persist for
several more days. The census of PV and N 2 shows that
symmetric instability must have played an active role
since there is a signiﬁcant volume of water that has zero
PV, but positive N 2.
There is no ambiguity about the depth of the mixed
layer when one considers the distribution of zonal-mean
potential vorticity (Fig. 13). One might choose the Q*
5 0.5 contour to distinguish between mixed layer and
interior ﬂuid, for example; the PV is normalized by the
initial value [Eq. (19b)]—the interior, undisturbed ﬂuid
has a value of 1. At day 3, before the baroclinic eddies
have evolved, there is good agreement between the one-
dimensional prediction and the mixed layer depth di-
agnosed by consideration of the PV ﬁeld. The PV in
this region is close to zero, as expected. By day 9, how-
ever, there is a signiﬁcant departure from the one-di-
mensional prediction of mixed layer depth. The zonal-
Unauthenticated | Downloaded 03/05/26 01:04 AM UTC


<\!-- Page 11 -->
644
VOLUME 28
J O U R N A L O F P H Y S I C A L O C E A N O G R A P H Y
FIG. 6. Vertical sections from the 2D integration (experiment 2) at day 9, for the central part of the channel. (a)
Temperature (contour spacing of 0.028C), (b) Ertel potential vorticity normalized by the PV of the initial condition, and
(c) tracer. In each ﬁgure the ﬂow and 1D mixed layer depth are shown as in Fig. 5. Peak speeds (y, w) are (0.11, 0.050)
m s21.
average PV shows that the deepening has been retarded
on the side of the channel where there is large surface
buoyancy loss, whereas it is increased on the weakly
cooled side. This systematic lateral ﬂux is provided by
baroclinic eddies and becomes the major contribution
to the buoyancy budget for columns of water in the
region where the cooling is weak. Comparing the angles
of mixed layer isotherms in the 2D and 3D integrations
Unauthenticated | Downloaded 03/05/26 01:04 AM UTC


<\!-- Page 12 -->
APRIL 1998
645
H A I N E A N D M A R S H A L L
FIG. 7. Two-dimensional switch-off experiment. As Fig. 6, but at day 10, after the surface cooling has ceased for 24
h. Absolute momentum surfaces are also shown by dotted lines (contour spacing 0.1 m s21). Peak speeds are (0.030,
0.0045) m s21 in the (y, z) directions at this time.
shows that the slopes are steeper in the 2D case (Figs.
6 and 11). This suggests that baroclinic instability is
more effective at transferring ﬂuid parcels horizontally
than symmetric instability. Below we study this hori-
zontal transfer process directly by considering the model
heat budget (Fig. 14).
The simulations shown in Figs. 8–14 exhibit consid-
erable nonhydrostatic behavior as measured by a nonhy-
Unauthenticated | Downloaded 03/05/26 01:04 AM UTC


<\!-- Page 13 -->
646
VOLUME 28
J O U R N A L O F P H Y S I C A L O C E A N O G R A P H Y
FIG. 8. The evolution of temperature at a depth of 65 m for days 3, 6, and 9 of experiment 3. The
switch over from small-scale gravitational convection to ﬁnite amplitude baroclinic instability and
geostrophic turbulence is clear.
drostatic parameter, n. Marshall et al. (1997a) show that
the condition for nonhydrostatic dynamics is that
n [ Ri Ro k 1,
where Ro is the Rossby number. We observe a transition
from nonhydrostatic plume dynamics through baroclinic
instability modiﬁed by nonhydrostatic effects to hydro-
static baroclinic instability as the lateral scale expands—
see Table 2 of Marshall et al. (1997a) where typical values
of n are presented in a range of numerical simulations that
span the hydrostatic to nonhydrostatic regime.
4. Buoyancy transfer by baroclinic eddies
Here we quantitatively assess the eddy transfer oc-
curing in our mixed layer and the way it modiﬁes the
Unauthenticated | Downloaded 03/05/26 01:04 AM UTC


<\!-- Page 14 -->
APRIL 1998
647
H A I N E A N D M A R S H A L L
FIG. 9. Sections of zonal mean temperature (full lines, contour spacing of 0.058C) and absolute mo-
mentum (dotted lines, contour spacing of 0.25 m s21) for days 3, 6, and 9 of experiment 3. On average
the mixed layer is weakly stratiﬁed by day 9 despite the ongoing surface buoyancy loss. The thick dashed
line in each plot is the prediction of the 1D law for the depth of the mixed layer [Eq. (5)]. Notice the x-
scale change from Figs. 5, 6, and 7.
development of the layer. Consider Fig. 14 showing the
evolution of the total heat content for the southern half
of the channel. The simple 1D prediction (surface heat
loss only) is plotted with the true model heat content
as they develop with time. The southern half of the
channel loses heat at a greater rate than the simple 1D
model, a consequence of the systematic lateral eddy
ﬂux; the northern half gains this heat through eddy trans-
Unauthenticated | Downloaded 03/05/26 01:04 AM UTC


<\!-- Page 15 -->
648
VOLUME 28
J O U R N A L O F P H Y S I C A L O C E A N O G R A P H Y
FIG. 10. (a) Horizontal section of potential vorticity (normalized by initial PV value) at 671-
m depth at day 9 from experiment 3, (b) volumetric census of PV and N 2 classes at day 9
(logarithmic scale to base 10) showing that a signiﬁcant fraction of the ﬂuid has negative PV
and N 2. The total volume of water in the integration is 3 3 1012 m3 and was used to normalize
the volumes in (b).
fer and so loses heat at a lesser rate than the 1D pre-
diction. The lateral eddy buoyancy ﬂux is diagnosed as
the difference between the 1D model and the evolution
observed in the model.
We now study the eddy ﬂux in the context of a closure
based on the local mean buoyancy gradient following
Green (1970) and Stone (1972). An alternative approach
is to ﬁnd the ﬂux required to return the ﬂow to a mar-
ginally stable state (Stone 1978)—the familiar analog
for gravitational instability is convective adjustment to
a neutrally stable state. However, recent work by Pavan
and Held (1996) and Vallis (1988) suggests that a bar-
oclinic adjustment procedure is less successful than a
local transfer relationship in predicting baroclinic eddy
ﬂuxes (at least in their two-layer quasigeostrophic beta-
plane models). Therefore, we pursue a gradient closure
and express the lateral buoyancy ﬂux across the channel,
y9b9, due to baroclinic instability as
y9b9 5 2Kby,
(20)
where K, an eddy transfer coefﬁcient, is to be related
to mean-ﬂow quantities. The lateral buoyancy gradient
is by and the overbar indicates a mean quantity (aver-
aged along the channel).
According to ‘‘mixing length’’ theory, the transfer
coefﬁcient can be expressed in terms of the character-
istic velocity and length scales of the transfer process,
thus
K 5 y9l9 5 cey eddym.
(21)
Here m is a measure of the lateral transfer length scale,
y eddy is a measure of the typical eddy velocity, and ce
Unauthenticated | Downloaded 03/05/26 01:04 AM UTC


<\!-- Page 16 -->
APRIL 1998
649
H A I N E A N D M A R S H A L L
is an efﬁciency factor. The velocity scale will be de-
duced from the thermal wind in the mixed layer, while
the transfer length scale will be related to the inherent
scales in the problem—to the deformation radius and
the width of the baroclinic zone.
Eddy velocity scale. We suppose, not unreasonably,
that the eddy velocity is given by
2
M h
y
5 |u| 5 u h 5
,
(22)
eddy
z
f
the maximum speed of the thermal wind in the mixed
layer. Here h is the mixed layer depth at the channel
center, M 2 5 |]b/]y| is a measure of the lateral strati-
ﬁcation, and the thermal wind relation has been used.
Transfer length scale. What is the appropriate spatial
scale to characterize the transfer process (mixing
length)? Two choices are readily apparent in our channel
calculation: the scale of the eddies themselves and the
width of the front on which they grow (controlled in
our calculations by the scale of the forcing function).
The strong impression one obtains by inspection of the
evolving ﬁelds of our numerical experiments (such as
Figs. 8 and 10) is that, as the eddies mature, they expand
to match the scale of the baroclinic zone.
We provisionally suppose, then, that the transfer
length scale is proportional to the frontal width. That
is,
m 5 Lzone.
(23)
In our experiments
Lzone 5 2Lf,
(24)
where Lf is the width of the hyperbolic tangent function
used to force the model [Eq. (18)].
Thus, using (23), we write (21) as
K 5 ceLzoneu,
(25)
where u 5 y eddy is given by (22).
Below we seek support for the form (25) from our
explicit calculations and determine the efﬁciency factor
ce.
Transfer timescale. Using the characteristic velocity
scale and the lateral transfer length scale allows one to
form a timescale for the transfer process itself. The
transfer timescale t transfer is
t transfer 5 m/y eddy ø Lzone/u
(26)
using (22) and (23).
If the transfer space scale is controlled by the defor-
mation radius rather than Lzone, then the above scaling
laws would be modifed by powers of the Burger number,
Bu,
Bu 5
/
,
2
2
L L
r
zone
(27)
where the deformation radius is
Lr 5 Nthh/ f.
The choice of transfer scale is one of the principle
features distinguishing the baroclinic eddy parameteri-
zations of Green (1970), who supposed that the transfer
scale was Lzone, and Stone (1972), who supposed that it
was Lr, the radius of deformation. This distinction is of
little consequence when the Burger number is nearly
constant, but it will turn out to be important here. In
our experiments (see Table 1) Bu varies by an order of
magnitude since Lr and Lzone can be speciﬁed indepen-
dently by changing the initial stratiﬁcation and the scale
of the externally imposed cooling function. We can thus
study the dependence of our closure ideas on the Burger
number. We now test the scaling law against our explicit
numerical calculations focusing on two central issues;
the space and timescale of the transfer process and the
closure for the magnitude of the eddy buoyancy ﬂux.
a. Testing the transfer timescale
Figure 14 shows the development of the cumulative
eddy heat ﬂux in one particular experiment—as time
goes by the eddies become more and more important.
We estimate from the model the time t model at which
lateral transfer by eddies has become signiﬁcant. We
then compare this to t transfer, Eq. (26), for each of the
14 experiments shown in Table 1.
The ﬁrst step is to relate t transfer to external parameters,
which are controlled in the experiments. Using (22), Eq.
(26) can be written:
L
f
zone
t
5
.
(28)
transfer
2
hM
To express M 2 in terms of external parameters we
assume that
1) the mixed layer smoothly connects with the ther-
mocline below,
M 2 5 |hy|
,
2
N th
(29)
and
2) the 1D law of mixed layer deepening applies at the
channel center [Eq. (5)].
Hence,
t
2
M 5 |B |N
.
(30)
y
th!2B
Substituting for M 2 in (28) using (30) and solving for
the time we ﬁnd
B f
t
5 2
.
(31)
transfer
2
!By
For the form of cooling chosen to drive our numerical
experiments, Eq. (18), the above expression can be eval-
uated at midchannel and expressed as
Lzone
t
f }
,
(32)
transfer
lrot
where lrot 5
B/ f 3 5 Ro*h is the rotational length
Ï
scale formed from the external parameters B and f eval-
Unauthenticated | Downloaded 03/05/26 01:04 AM UTC


<\!-- Page 17 -->
650
VOLUME 28
J O U R N A L O F P H Y S I C A L O C E A N O G R A P H Y
FIG. 11. Vertical sections from experiment 3 at day 9, half way along the channel (x 5 25 km). (a) Temperature
(contour spacing 0.028C), (b) Ertel potential vorticity normalized by the PV of the initial condition, and (c) tracer.
Velocities are as in Fig. 5 [peak speeds are (0.21, 0.022) m s21 in the (y, z) directions].
uated at midchannel; this scale, the physical interpre-
tation of which is discussed in Marshall et al. (1994),
crops up in many problems in rotating convection.
The instability timescale observed in the model mea-
sured against f, t model f, is plotted in Fig. 15 against Lzone/
lrot, Eq. (32), for 14 experiments in which the external
forcing, rotation rate, initial stratiﬁcation, and domain
size were all varied (Table 1). We see that, indeed, the
scaling law for t transfer f is a good prediction for the
timescale deduced from the model t model f. The numer-
ical results also conﬁrm that the transfer timescale is
independent of the ambient stratiﬁcation and the domain
Unauthenticated | Downloaded 03/05/26 01:04 AM UTC


<\!-- Page 18 -->
APRIL 1998
651
H A I N E A N D M A R S H A L L
FIG. 12. Three-dimensional switch-off experiment. As Fig. 10 but at day 10, 24 h after the
surface cooling has ceased. Both PV and N 2 are now nonnegative.
size. For example, experiments 6 and 8 were performed
with initial stratiﬁcations that varied by a factor of 4,
yet showed only an 11% change in the transfer times-
cale; experiments 3 and 16 were identical except that
experiment 16 was carried out in a channel twice as
wide, and the difference in transfer timescales was only
2% (see Table 1).
If, instead, we had chosen a transfer scale of Lr then,
replacing Lzone by Lr in Eq. (28), we ﬁnd that the transfer
timescale varies as (Lzone/Lrot)2/3. The implied curve is
also shown in Fig. 15 and is clearly less successful at
explaining the data than Eq. (32) assuming Lzone is the
transfer scale. It appears, then, that the appropriate trans-
fer scale is set by the width of the baroclinic zone, rather
than the instability scale.
b. The magnitude of the lateral ﬂuxes—Finding the
constant ce
Next we consider the magnitude of the eddy buoyancy
ﬂux driven by baroclinic eddies observed in our suite
of numerical experiments (Table 1) and the ability of
the simple ideas outlined above to predict it. We focus
on the cumulative eddy transfer of heat as a function
of time, E(t), plotted in Fig. 14. We suppose that E(t)
is due entirely to eddies since the plume convection/
symmetric instability is suppressed by the eddies (Fig.
8). For each integration we estimate ce by ﬁnding the
value that gives the best ﬁt to the time evolution of the
model using Eq. (25) to specify the value of K, evaluated
at t 5 t model. For example, the dashed line in Fig. 14
shows the result for experiment 1. The magnitude of ce
from all 14 experiments was 0.0817 6 0.023. We also
allowed the K to vary with time, changing it in pro-
portion to the observed speed of the thermal wind in
the mixed layer as suggested by (25) according to (22).
There is some change—in this case ce 5 0.119 6
0.031—since the majority of the eddy ﬂux occurs when
the eddies are mature.
We discuss the implication of our results for the pa-
rameterization of eddy ﬂuxes (and place them in the
context of other work on this subject) in section 5.
Unauthenticated | Downloaded 03/05/26 01:04 AM UTC


<\!-- Page 19 -->
652
VOLUME 28
J O U R N A L O F P H Y S I C A L O C E A N O G R A P H Y
FIG. 13. Sections of zonal mean Ertel potential vorticity (normalized by initial PV value) at days 3, 6, and 9 of
experiment 3. The thick dashed line in each plot is the prediction of the one-dimensional law for the depth of the mixed
layer [Eq. (5)].
5. Discussion
a. Implications for the oceanic mixed layer
How do our results apply to the mixed layer in a
general context? Here we address the issues of how they
may be modiﬁed in the presence of other processes, and
estimate their quantitative signiﬁcance.
Theory (section 2) suggests that a mixed layer with
spatial density gradients ought to convect along sloping
paths. This symmetric instability rapidly generates a lay-
er with vanishing potential vorticity (Ri 5 1), but non-
Unauthenticated | Downloaded 03/05/26 01:04 AM UTC


<\!-- Page 20 -->
APRIL 1998
653
H A I N E A N D M A R S H A L L
FIG. 14. Evolution of terms in the heat budget for the southern half
of the channel in experiment 3. The circles show the evolution of
the observed heat loss—it increases with time due to the surface
cooling (shown by the dashed line) and eddy heat transfer to the
north. The full line shows the prediction of the eddy parameterization
scheme as explained in the text, where ce 5 0.0827.
FIG. 15. Regression of normalized transfer timescale measured
from the model versus Lzone/lrot. Each symbol is the result of a different
experiment from Table 1 (experiments 3–16). The full line has a
gradient of unity, corresponding to the scaling law given by Eq.
(25)—where the transfer scale is the width of the baroclinic zone.
The dashed line has a gradient of 2/3, corresponding to the case in
which the transfer scale is the deformation radius. The t model is
deﬁned to be the time at which the cumulative heat lost from the
southern half of the channel is twice the cumulative heat lost from
the surface (Fig. 14).
zero vertical stratiﬁcation. Thereafter, a nonhydrostatic
baroclinic instability develops that provides a lateral and
vertical buoyancy ﬂux that becomes the dominant mode
of buoyancy transfer. The fully nonlinear numerical ex-
periments are in general support, and indicate that the
PV of the surface-forced layer is reset to zero on a
timescale of a few hours, with pronounced stratiﬁcation
in the overturning ﬂuid. Baroclinic eddies subsequently
develop, and we suggest scaling laws for the magnitude
of the buoyancy transfer, that found support in the nu-
merical model.
The model does not, of course, attempt to include all
the processes at work in the upper ocean, however. There
are processes other than spatial gradients in the buoy-
ancy forcing or stratiﬁcation that cause lateral density
gradients in the mixed layer. The real ocean develops
surface frontal regions from larger-scale strain (due to
instability of the main thermocline for example) and
outcropping of density surfaces. Lateral contrasts in me-
chanical forcing that result from isolated atmospheric
disturbances and remotely produced waves can also be
important. These include mixing due to breaking surface
waves generated by the wind, inertial waves, and Lang-
muir turbulence. Upwelling at a coastal boundary is
another important source of nearshore fronts.
For the purpose of simplicity, and to facilitate un-
derstanding, we have not attempted to represent these
processes in our model. They may modify or disrupt
the thermal wind that is ultimately set up in the mixed
layer. But it is reasonable to suppose that deep mixed
layers in which lateral gradients persist for several ro-
tation periods will come close to a balanced state. In-
deed, one of the key ﬁndings of our study is that currents
in thermal wind balance in the mixed layer, averaged
over a few kilometers and days, do persist despite strong
surface forcing. The consequence is a vigorous baro-
clinic instability. Eddies develop over several days—
somewhat slower than the timescale for convection/
symmetric instability—and efﬁciently redistribute
buoyancy within the mixed layer.
What is the signiﬁcance of the mixed layer baroclinic
instability for the real ocean? We will use our param-
eterization derived in section 4 and calibrated against
our numerical calculations, to make inferences from ob-
servations of sea surface density about the magnitude
of the transfer coefﬁcient K. Combining (25) and (22),
we obtain
gL
|s |h
zone
y
K 5 c L
|b |h/ f 5 c
,
(33)
e
zone
y
e
r f
0
with ce set equal to 0.1.
It is unclear what Lzone should be, but let us assume
a conservatively low estimate of 20 km based on the
FASINEX survey (Pollard and Regier 1992). We can
now use estimates of sy from observations of sea sur-
face density made on a 10 000-km cruise track in the
northeast Atlantic during April and May 1991 shown
in Fig. 16a. These data were taken with a thermosali-
Unauthenticated | Downloaded 03/05/26 01:04 AM UTC


<\!-- Page 21 -->
654
VOLUME 28
J O U R N A L O F P H Y S I C A L O C E A N O G R A P H Y
FIG. 16. (a) Sea surface density (kg m23) measured every 250 m along a cruise track in the northeast
Atlantic in April and May 1991. The cruise involved six north–south legs between 428 and 548N from
the European continental slope to the Mid-Atlantic Ridge (Cunningham et al. 1992). Mixed layer depths
ranged from around 50 m to a few 100 m (Fig. 2 shows the upper oceanic structure from the ﬁrst 300
km of the cruise track). Seawater was drawn on board from a depth of 6 m and analyzed using a
thermosalinograph then averaged to give a sea surface density value every 250 m. (b) Distribution of sea
surface density gradients (kg m24). The data in (a) have been smoothed using a low-pass ﬁlter with
wavelength 5 km. (c) Distribution of transfer coefﬁcients (m 2 s21) calculated from Eq. (33) using the
density gradient distribution from (b). Four mixed layer depths h have been assumed (50, 100, 200, and
400 m). Also shown is a constant transfer coefﬁcient of 250 m2 s21 typical of the values used in oceanic
general circulation models.
Unauthenticated | Downloaded 03/05/26 01:04 AM UTC


<\!-- Page 22 -->
APRIL 1998
655
H A I N E A N D M A R S H A L L
nograph, giving an average measurement every 250 m
along the cruise track. Before computing the gradient,
the dataset was smoothed with a low-pass ﬁlter, cutting
off signals shorter than 5 km. The gradient was binned
and its distribution is shown in Fig. 16b. The mixed
layer depth varied between several tens of meters and
a few hundreds (for example, the data in Fig. 2, from
the ﬁrst part of the cruise, shows depths of 200–400
m), so we compare results for four characteristic depths
(50, 100, 200, and 400 m).
Figure 16c shows the estimate of transfer coefﬁcient
due to mixed layer baroclinic instability. We plot the
distribution of K corresponding to each range of sy for
the four mixed layer depth choices. Typically the trans-
fer coefﬁcients are in the range from 100 to 1000 m2
s21 and are largest for deep mixed layers and strong
fronts. [Typical values of eddy diffusivity used in eddy-
resolving basin models are O(100 m2 s21).] It is clear
that a constant K (constant Fickian diffusion—com-
monly used in oceanic general circulation models) poor-
ly represents the eddy transfer process; the coefﬁcient
estimated using (33) varies by three orders of magnitude
for a given depth. We have repeated these estimates
using a 1-km low-pass ﬁlter and ﬁnd that the transfer
coefﬁcients are increased by a factor of 10 in this case.
Our calculations here (and elsewhere—see, e.g.,
Jones and Marshall 1993; Visbeck et al. 1996) clearly
show that lateral transfer by mixed layer baroclinic ed-
dies is important on the margins of deep convection
chimneys. But our calculations also suggest that they
may be important in less extreme mixing regimes such
as frontal regions associated with O(100 m) deep mixed
layers. It is reasonable to speculate that baroclinic in-
stability of the mixed layer is commonplace and pro-
vides a signiﬁcant and efﬁcient mechanism of lateral
and vertical buoyancy transport through them. The sea-
sonal cycle modulates the process leading to greatest
power when there is deepest mixing. As such, it will
be most signiﬁcant when the ocean interior receives
water fresh from atmospheric contact and so will be
responsible, in part, for determining the character of
newly formed thermocline water.
b. Implications for eddy parameterizations
Green (1970) and Larichev and Held (1996) suggest
the following closure for the eddy transfer coefﬁcent:
2
m
K
5 c
,
(34)
G
Gt Eady
where t Eady 5
Ri/ f 5 Lr/u is a measure of the in-
Ï
stability timescale [Eq. (15)] and m is the transfer scale.
Green chose m 5 Lzone, the width of the baroclinic zone.
Instead, Stone (1972) supposed that
KS 5 cSum
(35)
as in (25), arguing that u should be given by (22) and
setting m 5 Lr, the deformation radius.
The forms (34) and (35) should be compared to K 5
ceLzoneu, Eq. (25), the form that best ﬁts the data here.
It is straightforward to show
ce 5 cS
Bu 5 cG/
Bu,
Ï
Ï
where cS and cG are the empirical constants used in the
Stone and Green theories. The average Burger number
for the experiments performed here is 0.063, so the ce
ø 0.1 found in the present study implies cS 5 0.42 and
cG 5 0.025. These are in acceptable agreement with
Stone’s (1972) estimate of 0.86 and the Visbeck et al.
(1996) value for cG of 0.015. Clearly, one can only
distinguish between these closure theories by varying
the Burger number, as we have done here. We ﬁnd that
(25), the simplest of the three, gives the best ﬁt to the
data.
In contrast to the problem we have addressed in this
paper, several recent studies on parameterizing baro-
clinic eddy ﬂuxes have focused on two-layer quasi-
geostrophic models (Pavan and Held 1996; Larichev and
Held 1995; Held and Larichev 1996; Vallis 1988). Lar-
ichev and Held (1995) ﬁnd support for Green’s scaling
in statistically steady models of homogeneous turbu-
lence, while Held and Larichev (1996) conﬁrm that the
closure still applies on a beta plane, provided the Rhines
(1975) length scale is used for m in Eq. (34). Our model
differs from these studies in two important respects:
First, we solve the nonhydrostatic Navier–Stokes equa-
tions for ﬂows that are not quasigeostrophic (clearly
evident from the large vertical isotherm excursions in
Fig. 11). Second, there is no statistical steady state and
our eddy ﬁeld is inhomogeneous due to the shape of
the imposed forcing function. Despite these limitations
we are encouraged by the success of the scaling we
propose. However, we consider it provisional and it must
be tested against statistically steady, homogeneous sim-
ulations of baroclinic instability.
Acknowledgments. We thank Martin Visbeck for his
helpful and enthusiastic remarks and the comments of
two anonymous referees. The datasets were taken from
the United Kingdom’s Vivaldi ’91 cruise. TWNH and
JCM were supported by grants from the ACCP program
of NOAA (NA46GP0125) and Grant (OCE-9503895)
of NSF.
APPENDIX A
Parcel Theory of Gravitational Instability
Suppose that a weakly stratiﬁed ocean is subject to
vigorous cooling at the surface over some hundreds of
kilometers, producing a density inversion and the pos-
sibility of overturning. The ﬂuid cannot simultaneously
overturn on this scale; rather, the qualitative response
to widespread cooling is one in which relatively small
convection cells (plumes) develop. The detailed physics
setting the plume scale is as yet unclear, but clearly the
Unauthenticated | Downloaded 03/05/26 01:04 AM UTC


<\!-- Page 23 -->
656
VOLUME 28
J O U R N A L O F P H Y S I C A L O C E A N O G R A P H Y
gross transfer properties of the population of convective
cells must be controlled by the large scale; the vertical
buoyancy ﬂux due to convection must offset loss at the
surface. A law of vertical buoyancy transfer for the
plume scale and its instability timescale can be devel-
oped using parcel theory as follows.
Suppose that the net effect of overturning is to ex-
change particles of ﬂuid, of density r1 and r2, over a
depth Dz. Dense water sinks displacing lighter water
below, and releases potential energy to power the plume
and ﬂux buoyancy vertically. The change in potential
energy DP consequent on this idealized rearrangement
of particles (Fig. 1a) is given by
DP 5 P
2 P
final
initial
5 g[(r z 1 r z ) 2 (r z 1 r z )]
1 2
2 1
1 1
2 2
5 r DbDz,
(A1)
0
where Db 5 (b2 2 b1) is the buoyancy difference of
the particles exchanged over a distance Dz 5 (z2 2 z1),
g is the acceleration due to gravity, and r0 is a repre-
sentative value of the density.
Equating the released potential energy to the acquired
kinetic energy of the ensuing convective motion, DK 5
2 3 (1/2)r0w2 (where w is the vertical velocity scale
and there is a factor of 2 because there are two particles),
then
w2 5 DbDz.
(A2)
Now, if heavy ﬂuid lies above light ﬂuid, an unstable
disturbance with growth rate v will develop. Thus, set-
ting z } evt; d/dt 5 v; w 5 dz/dt 5 vz, then (A2)
implies
2
2
2
2
v z 5 |N |z ,
so
2
2
v 5 |N |,
(A3)
which is the result given by linear stability analysis [in
the inviscid limit, the lateral scale of the fastest growing
convective mode collapses to zero, no energy is in-
volved in lateral motion, and the limit (A3) is achieved].
The implied vertical buoyancy ﬂux on the plume scale
is then, using (A2),
Bp 5 wDb 5 Dz1/2Db3/2,
(A4)
which is Eq. (4).
APPENDIX B
Energy Analysis of the Thermal Wind
We now consider parcels of incompressible ﬂuid at
positions (y1, z1) and (y2, z2) in thermal wind balance
with a meridional density gradient. The parcels are then
interchanged adiabatically (i.e., with conservation of
buoyancy).
a. Potential energy
The change in potential energy (per unit volume) is
again (A1):
DP 5 r (z 2 z )(b 2 b ).
(B1)
o
2
1
2
1
We also have
]b
]b
Db 5
Dy 1
Dz
]y
]z
so that
2
2
(b 2 b ) 5 M (y 2 y ) 1 N (z 2 z )
2
1
2
1
2
1
if, for simplicity, N 2 and M 2 are assumed constant and
measure the strength of the vertical and horizontal den-
sity gradients, respectively.
The slope sb of a buoyancy surface is, setting b2 5
b1,
2
dy
M
s 5
5 2
.
b
2
dz
N
Hence, if s 5 (z2 2 z1)/(y2 2 y1) is the slope of the
surface along which the particles are interchanged, we
may write
b2 2 b1 5 N 2(y2 2 y1)(s 2 sb),
so (B1) becomes
DP 5 r0N 2Dy2s(s 2 sb).
(B2)
In convectively stable conditions N 2 . 0. Hence, the
sign of DP is the same as that of the factor s(s 2 sb)
and it will be negative, corresponding to the possibility
of instability, if s , sb, that is, if the slope of the ex-
change surface has a smaller slope than that of the iso-
pycnals.
b. Kinetic energy
We must now consider the change in kinetic energy
involved in exchanging ﬂuid parcels. Let us assume a
zonal motion u 5 u(y, z) of zonal tubes of ﬂuid inde-
pendent of x (rather than just parcels) so that zonal
momentum cannot be changed by pressure gradient
forces. Then we have
Dm 5 0,
Dt
where m 5 u 2 fy is the absolute momentum. For any
small displacement in the y direction we must have
Du 5 fDy
as the change following the motion.
Now consider the change in kinetic energy resulting
from the exchange of tubes of ﬂuid with zonal motion
u1 at (y1, z1) and u2 at (y2, z2). Then
Unauthenticated | Downloaded 03/05/26 01:04 AM UTC


<\!-- Page 24 -->
APRIL 1998
657
H A I N E A N D M A R S H A L L
1
2
DK 5
r [{u 1 f (y 2 y )}
0
1
2
1
2
2
2
2
1 {u 2 f (y 2 y )} 2 u 2 u ]
2
2
1
1
2
(u 2 u )
2
1
2
5 r (y 2 y ) f f 2
.
0
2
1 [
]
(y 2 y )
2
1
Since
]u
]u
(u 2 u ) 5
(y 2 y ) 1
(z 2 z ),
2
1
2
1
2
1
]y
]z
this may be written
]u
]u
2
DK 5 r (y 2 y ) f
f 2
2 s
,
0
2
1 1
2
]y
]z
or, since f ]u/]z 5 2M 2 5 N 2sb, then
]u
2
2
DK 5 r Dy
f
f 2
2 N ss .
(B3)
0
b
1
2
[
]
]y
Hence the change in total energy of the mean motion,
DE, is, adding (B3) to (B2) and rearranging,
2


]u
2
f


1 2
]z


]u
2
2
2
DE 5 r Dy
f
f 2
2
1 N (s 2 s ) .
0
b


2
1
2
]y
N


c. Symmetric instability
Regarding s, the direction of exchange, as a variable,
DE has a minimum when s 5 sb, that is, when the
exchange is in the initial isopycnal surface. Then
z
1
2
2
(DE)
5 r Dy f
2
(B4)
min
0
1
2
f
Ri
will be negative if
f
Ri ,
,
(B5)
z
corresponding to negative Ertel potential vorticity and
the possibility of instability. Here z 5 f 2 ]u/]y is the
vertical component of absolute vorticity and
2
2
f N
Ri 5
(B6)
4
M
is the Richardson number.
From (B4) we may again deduce the growth rate:
z
1
2
2
v 5 2 f
2
.
(B7)
1
2
f
Ri
d. Baroclinic instability
If the Richardson number of the mean ﬂow is large,
then only changes in its potential energy (B2) need be
taken into account. For a given exchange distance, (y2
2 y1), the release of potential energy will be a maximum
(2DP a maximum) when s(s 2 sb) is a maximum; in
other words, when
sb
s 5
,
2
then
1
2
2
2
(DP)
5 2 r N Dy s
min
0
b
4
2
1
f
2
5 2 r Dy
.
0
4
Ri
We see that in any region where a thermal wind exists
it is always possible to release potential energy for eddy
growth provided an appropriate rearrangement of par-
ticles takes place.
If an unstable disturbance grows, then y } evt; d/dt
5 v; y 5 dy/dt 5 vy. Equating released potential energy
to acquired kinetic energy of the eddying motion (r0y 2)
we ﬁnd, in direct analogy with the upright convection
problem outlined above,
2f
2
v ø
,
(B8)
Ri
which is a heuristic derivation of the growth rate of an
Eady wave. Equation (B8) is in effect just (A3) but with
the |N 2| measured along a slope that has one-half that
of the isopycnals. Parcel theory has led us to the same
result one derives from linear stability analysis.
REFERENCES
Barth, J. A., 1994: Short-wavelength instabilities on coastal jets and
fronts. J. Geophys. Res., 99, 16 095–16 115.
Chandrasekar, S., 1961: Hydrodynamic and Hydromagnetic Stability.
Dover, 652 pp.
Cunningham, S. A., and Coauthors, 1992: SeaSoar CTD, ﬂuorescence
and scalar irradiance data from RRS Charles Darwin cruises
58/59, NE Atlantic (Vivaldi 91). IOSDL Rep. 299, 48 pp. [Avail-
able from Library, SOC, Empress Dock, Southampton SO14
3ZH, United Kingdom.]
Eady, E. T., 1949: Long waves and cyclone waves. Tellus, 1, 33–52.
, 1951: The quantitative theory of cyclone development. Com-
pendium of Meteorology, T. F. Malone, Ed., Amer. Meteor. Soc.,
464–469.
Emanuel, K., 1994: Atmospheric Convection. Oxford University
Press, 580 pp.
Fukamachi, Y., J. P. McCreary, and J. A. Proehl, 1995: Instability of
density fronts in layer and continuously stratiﬁed models. J.
Geophys. Res., 100, 2559–2577.
Green, J. A., 1970: Transfer properties of the large-scale eddies and
the general circulation of the atmosphere. Quart. J. Roy. Meteor.
Soc., 96, 157–185.
Held, I. M., and V. D. Larichev, 1996: A scaling theory for hori-
zontally homogeneous baroclinically unstable ﬂow on a beta
plane. J. Atmos. Sci., 53, 946–952.
Unauthenticated | Downloaded 03/05/26 01:04 AM UTC


<\!-- Page 25 -->
658
VOLUME 28
J O U R N A L O F P H Y S I C A L O C E A N O G R A P H Y
Hoskins, B. J., 1974: The role of potential vorticity in symmetric
stability and instability. Quart. J. Roy. Meteor. Soc., 100, 480–
482.
Jones, H., and J. C. Marshall, 1993: Convection with rotation in a
neutral ocean: A study of open-ocean deep convection. J. Phys.
Oceanogr., 23, 1009–1039.
Kraus, E. B., and J. S. Turner, 1967: A one-dimensional model of
the seasonal thermocline. II The general theory and its conse-
quences. Tellus, 19, 98–105.
Larichev, V. D., and I. M. Held, 1995: Eddy amplitudes and ﬂuxes
in a homogeneous model of fully developed baroclinic insta-
bility. J. Phys. Oceanogr., 25, 2285–2297.
Lin, S.-J., and R. T. Pierrehumbert, 1988: Does Ekman friction sup-
press baroclinic instability? J. Atmos. Sci., 45, 2920–2933.
Marshall, J. C., J. Whitehead, and T. Yates, 1994: Laboratory and
numerical experiments in oceanic convection. Ocean Processes
in Climate Dynamics, P. Malanotte-Rizzoli and A. Robinson,
Eds., Kluwer Academic Press, 437 pp.
, C. N. Hill, L. Perelman, and A. Adcroft, 1997a: Hydrostatic,
quasi-hydrostatic and non-hydrostatic ocean modeling. J. Geo-
phys. Res., 102, 5733–5752.
, A. Adcroft, C. N. Hill, L. Perelman, and C. Heisey, 1997b: A
ﬁnite-volume, incompressible Navier Stokes model for studies
of the ocean on parallel computers. J. Geophys. Res., 102, 5753–
5766.
Mellor, G. L., and T. Yamada, 1974: A hierarchy of turbulence closure
models for planetary boundary layers. J. Atmos. Sci., 31, 1791–
1806.
Pavan, V., and I. M. Held, 1996: The diffusive approximation for
eddy ﬂuxes in baroclinically unstable jets. J. Atmos. Sci., 53,
1262–1272.
Pollard, R. T., 1986: Frontal surveys with a towed proﬁling conduc-
tivity/temperature/depth measurement package (SeaSoar). Na-
ture, 323, 433–435.
, and L. A. Regier, 1992: Vorticity and vertical circulation at an
ocean front. J. Phys. Oceanogr., 22, 609–625.
, P. B. Rhines, and R. O. R. Y. Thompson, 1973: The deepening
of the wind-mixed layer. Geophys. Fluid. Dyn., 3, 381–404.
Price, J. F., R. A. Weller, and R. Pinkel, 1986: Diurnal cycling: Ob-
servations and models of the upper ocean response to diurnal
heating, cooling and wind mixing. J. Geophys. Res., 91, 8411–
8427.
Rayleigh, O. M., 1916: On convection currents in a horizontal layer
of ﬂuid, when the higher temperature is on the lower side. Philos.
Mag. Ser. 6, 32, 529–546.
Rhines, P. B., 1975: Waves and turbulence on a beta-plane. J. Fluid
Mech., 69, 417–443.
Samelson, R. M., and C. A. Paulson, 1988: Towed thermistor chain
observations of fronts in the subtropical north Paciﬁc. J. Geo-
phys. Res., 93, 2237–2246.
Stone, P. H., 1971: Baroclinic instability under non-hydrostatic con-
ditions. J. Fluid Mech., 45, 659–671.
, 1972: A simpliﬁed radiative–dynamical model for the static
stability of rotating atmospheres. J. Atmos., Sci., 29, 405–418.
, 1978: Baroclinic adjustment. J. Atmos., Sci., 35, 561–571.
Turner J. S., 1973: Buoyancy Effects in Fluids. Cambridge University
Press, 368 pp.
Vallis, G. K., 1988: Numerical studies of eddy transport properties
in eddy-resolving and parametrized models. Quart. J. Roy. Me-
teor. Soc., 114, 183–204.
Veronis G., 1958: Cellular convection with ﬁnite-amplitude in a ro-
tating ﬂuid. J. Fluid Mech., 5, 410–435.
Visbeck, M., J. C. Marshall, and H. Jones, 1996: On the dynamics
of convective ‘‘chimneys’’ in the ocean. J. Phys. Oceanogr., 26,
1721–1734.
Young, W. R., and L. Chen, 1995: Baroclinic instability and ther-
mohaline gradient alignment in the mixed layer. J. Phys. Ocean-
ogr., 25, 3172–3185.
Unauthenticated | Downloaded 03/05/26 01:04 AM UTC
