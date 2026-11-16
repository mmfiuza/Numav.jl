# Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

using Luxor

# Julia-style palette, plus blue to complete the tetrahedron
const BLUE   = (0.251, 0.388, 0.847)  # #4063D8
const RED    = (0.796, 0.235, 0.200)  # #CB3C33
const GREEN  = (0.220, 0.596, 0.149)  # #389826
const PURPLE = (0.584, 0.345, 0.698)  # #9558B2

filename = "docs/src/assets/logo"
D = 100
L = D * sqrt(2)/2
margin = 2
    
# Vertex positions (front triangle + back vertex tucked toward center)
p1 = Point(0, 0)
p2 = Point(L/2, L*sqrt(3)/2)
p3 = Point(D-L, D)
p4 = Point(L, 0)
p5 = Point(D, D)
p6 = Point(D-L/2, D-L*sqrt(3)/2)

Drawing(D + margin, D + margin, "docs/src/assets/logo.svg")
translate(margin/2, D+margin/2)
scale(1, -1)

# Faces
sethue(RED...)
poly([p1, p2, p3], :fill, close=true)
sethue(PURPLE...)
poly([p2, p3, p6], :fill, close=true)
sethue(GREEN...)
poly([p2, p4, p6], :fill, close=true)
sethue(BLUE...)
poly([p4, p5, p6], :fill, close=true)

# Edges
setdash("solid")
setlinecap("round")
setline(1)
sethue(0, 0, 0)
line(p1, p2, :stroke)
line(p2, p3, :stroke)
line(p2, p4, :stroke)
line(p3, p1, :stroke)
line(p2, p6, :stroke)
line(p3, p6, :stroke)
line(p4, p5, :stroke)
line(p4, p6, :stroke)
line(p5, p6, :stroke)

finish()
