#!/usr/bin/python
def box(com_2): # com is list, xyz is int
	
	x = 1
	y = 1
	z = com_2[2]
	xs = []
	ys = []
	zs = []
	points = []
	p1 = [com_2[0]+x, com_2[1], com_2[2]+15]
	points.append(p1)
	p2 = [com_2[0]-x, com_2[1], com_2[2]+15]
	points.append(p2)
	p3 = [com_2[0], com_2[1]+y, com_2[2]+15]
	points.append(p3)
	p4 = [com_2[0], com_2[1]-y, com_2[2]+15]
	points.append(p4)
	p5 = [com_2[0]+x, com_2[1], com_2[2]-15]
	points.append(p5)
	p6 = [com_2[0]-x, com_2[1], com_2[2]-15]
	points.append(p6)
	p7 = [com_2[0], com_2[1]+y, com_2[2]-15]
	points.append(p7)
	p8 = [com_2[0], com_2[1]-y, com_2[2]-15]
	points.append(p8)
	
	for point in points:
		
		xs.append(point[0])
		ys.append(point[1])
		zs.append(point[2])
		
	return xs,ys,zs

