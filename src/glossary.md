# Glossary

- entity: Some geometrical body (points, lines, surfaces, volumes, etc).

# Acronyms:

- IDX (index): An unsigned integer that uniquely identifies an entity.

- er (external region): An unsigned integer present in an external mesh file (.bdf, etc) or API that classifies an entity into some group, usually a material. An er sequence is NOT necessarily consecutive, since mesh formats can have arbitrary conventions.
    * esr: external surface region
    * evr: external volume region

- ir (internal region): An unsigned integer that classifies an entity into some group in Numav's internal implementation. An ir sequence is ALWAYS consecutive starting at 0, because that is conveninent to manipulate in C++.
    - isr: internal surface region
        - isri: internal surface region with impedance
        - isrv: internal surface region with velocity
        - isrp: internal surface region with pressure
    - ivr: internal volume region

- fi (frequency independent): Some quantity that will be the constant relative to sound frequency.

- fd (frequency dependent): Some quantity that will change relative to sound frequency.

- nnz (non-zero): Some quantity that is not null. Usually refers to the non-zero terms of a sparse matrix.

- sei (surface element index):
    - ISEI (impedance surface element index)
    - VSEI (velocity surface element index)
    - PSEI (pressure surface element index)
- vei (volume element index)

- fipi (frequency independent part index)

- gpi (gauss point index)

- di (dimension index)

- ni (node index)

- nci (node combination index)

- eni (elemental node index)
    - enis (elemental node index of surface)
    - eniv (elemental node index of volume)

- vpi (velocity point index)

- ppi (pressure point index)

- pni (pressure node index)

- pvi (pressure value index)

- nzi (non-zero index)

- pq (physical quantity)

- pqv (physical quantity value)
