README
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

<img src='man/figures/logo.png' align="right" height="139" />

# `lbmech`: <br/>The Mechanics of <br/>Landscape and Behavior

<!-- badges: start -->
<!-- badges: end -->

`lbmech` is a geospatial package intended to study the mechanics of
landscape and behavior. The entire project has been conceived of in
seven parts, of which only four are in mature-enough states for public
presentation/release. As such, the package should be considered largely
a lengthy work-in-progress. The functions are very heavily documented
under the ‘References’ section—each function contains at least one
minimum working example needed to generate usable data for every entry.

The chapters under the ‘Articles’ section provide extensive theoretical
grounding behind the functions in each relevant part, before
demonstrating their use in a non-trivial case usually involving some
form of large data analysis. They can be considered as a series ‘living’
working papers; the aim is to eventually publish each of them
independently, yet retain them together as part of an ever-growing
corpus of analytical approaches. Use of the code or functions discussed
in each part should use the appropriate citation.

Backwards-compatibility with previous versions cannot be guaranteed
until the release of version 1.0 unless explicitly stated; once
particular parts have been accepted post peer review, they will be
labeled as `completed`’ and the corresponding functions will be
guaranteed to continue to behave as intended for future versions.

## Part 1: Movement

[Part 1 is a set of highly-efficient functions for GIS-style
cost-distance
analysis.](https://andresgmejiar.github.io/lbmech/articles/movement.html)
It is currently in `gamma`, and therefore backwards compatibility cannot
yet be guaranteed, although the functions are robust and significantly
more efficient than
[`gdistance`](https://agrdatasci.github.io/gdistance/) and stability is
expected. The syntax is expected to remain consistent enough that
plotting functions have been developed and are included.

It was originally designed to allow for the calculation of time and
energetic/thermodynamic costs when moving across the landscape. The
default parameters allow for such calculations, but currently any
arbitrary cost function depending on a single raster input is supported
(Previous versions \< 0.4.0 supported any arbitrary cost function with
any number of raster inputs; this functionality will be returned in a
coming release). Version 0.2.0 was introduced in Mejia Ramon (2021:
136-240, 300-327) and was most recently presented at Mejia Ramon (2022).

## Part 2: Inequality

[Part 2 is a set of inferential tools to detect spatial inequality in
vector
datasets.](https://andresgmejiar.github.io/lbmech/articles/inequality.html).
It is currently in `gamma`, and therefore backwards compatibility cannot
yet be guaranteed, although the functions are quite efficient and
stability is expected. Moreover, they are the only robust method for
inference regarding various types of within-and between group inequality
in spatial contexts. The syntax is expected to remain consistent enough
that plotting functions have been developed and are included.

Its current implementation is limited to the set of error-based
inequality metrics (such as the Gini and Inoua indexes), however the
approach is applicable to a number of other spatial statistics including
for autocorrelation. Version 0.2.0 first presented in Munson et
al. (2023).

## Part 3: Productivity

Part 3 is a set of tools to interpolate agricultural productivity data.
It is currently in early `beta`. Therefore, while the functions are
included and heavily documented, a theory and applications chapter is
not yet included nor have they been rigorously stress tested beyond
their intended applications. Moreover, the base data has not yet been
included in the compiled package but will be shortly.

The current implementation is designed to take municipal-level
aggregates of agricultural productivity for a number of crops in the
various agricultural censuses of mid-1900s Mexico and convert them into
location-level predictions of expected productive quantiles. Version
0.2.0 was introduced in Mejia Ramon (2021: 71-135, 284-299) and was most
recently presented at Mejia Ramon et al. (2023).

## References

Inoua, Sabiou (2021). “Beware the Gini Index! A New Inequality Measure.”
*ESI Working Paper* 21-18.
<https://digitalcommons.chapman.edu/esi_working_papers/355/>

Mejia Ramon, Andres G. (2021). *Agricultural Productivity and
Human-Landscape Dynamics in the Early Basin of Mexico*.
Ph.D. dissertation in Anthropology, The Pennsylvania State University,
University Park. <https://doi.org/10.26207/rj25-0380>

Mejia Ramon, Andres G.(2022). “The Mechanics of Landscape and Behavior:
A Least-Cost Path Approach Grounded in Mechanical Physics.” Paper
presented at the *87th Annual Meeting of the Society for American
Archaeology*, Chicago.

Mejia Ramon, Andres G., Jessica L. Munson, Jill Onken, and Lorena Paiz
Aragon (2023). “Let the Crops Speak for Themselves: How to Avoid
Imposing Agroecological Assumptions at Altar de Sacrificios.” Paper
presented at the *87th Annual Meeting of the Society for American
Archaeology*, Chicago.

Munson, Jessica L., Andres G. Mejia Ramon, Lorena Paiz Aragon, Jill
Onken, and Jonathan Scholnick (2023). “Settlement Density, Household
Inequality, and Social Interaction in the Western Maya Lowlands.” Paper
presented at the *88th Annual Meeting of the Society for American
Archaeology*, Portland, Oregon.
