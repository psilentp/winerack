__author__ = 'psilentp'

# -*- coding: utf-8 -*-
# <nbformat>2</nbformat>

import svgwrite
from svgwrite import mm
import numpy as np
from numpy import pi,cos,sin,tan,dot
from numpy import pi,cos,sin,tan,dot,degrees
laser_kerf = 0.15
inside_rad = 55.0
wall_thick = 20.0 + (2 * laser_kerf)
cutout_depth = wall_thick/2.0 - laser_kerf
material_thick = 6.95 ###check me
sep = (304.8 - material_thick)/3.0
front_thick = 2*material_thick
dowel_rad = 6.45/2 - laser_kerf#3/8"
hex_y = lambda rad: rad*tan(pi/6)
hex_rad = lambda rad: np.linalg.norm(np.array([rad,hex_y(rad)]))
wall_len = lambda rad: 2*rad*sin(pi/6)
shift_dist = lambda d: d / (1+ 1/(cos(pi/6.0)))
shift = shift_dist(hex_rad(inside_rad+wall_thick) - hex_rad(inside_rad))*cos(pi/6)
num_rows = 2
num_cols = 4
num_nubins = 2
nubin_d = 0.10


def to_hex_coord(row,col):
    row_mod = (row%2)*(inside_rad+wall_thick)
    xv = np.array([((inside_rad+wall_thick))*2,0])
    yv = np.array([0,(inside_rad+wall_thick)*2*cos(pi/6)])
    return (np.array([col,row])*(xv+yv)) + np.array([row_mod,0])


def make_circles(p_rad,c_rad):
    p1 = np.array([p_rad,hex_y(p_rad),0])
    r_mat = rotation_matrix(np.array([0,0,1]),pi/3)
    p2 = dot(r_mat,p1)
    p3 = dot(r_mat,p2)
    p4 = dot(r_mat,p3)
    p5 = dot(r_mat,p4)
    p6 = dot(r_mat,p5)
    pvecs = [(point[0],point[1]) for point in [p1,p2,p3,p4,p5,p6]]
    circles = list()
    for p in pvecs:
        #tran = totrans("x+%s*2,y+%s"%(p[0],p[1]))
        circles.append(dwg.circle(fill='none',
            stroke='black',
            stroke_width = 1,
            r=c_rad,center = (p[0],p[1]) ))
    return circles

def make_hexes(p_rad,h_rad):
    p1 = np.array([p_rad,hex_y(p_rad),0])
    r_mat = rotation_matrix(np.array([0,0,1]),pi/3)
    p2 = dot(r_mat,p1)
    p3 = dot(r_mat,p2)
    p4 = dot(r_mat,p3)
    p5 = dot(r_mat,p4)
    p6 = dot(r_mat,p5)
    pvecs = [(point[0],point[1]) for point in [p1,p2,p3,p4,p5,p6]]
    hexes = list()
    for p in pvecs:
        #tran = totrans("x+%s*2,y+%s"%(p[0],p[1]))
        h = make_hex(h_rad)
        h.translate(tx = p[0],ty = p[1])
        hexes.append(h)
    return hexes


def make_hex(rad):
    p1 = np.array([rad,hex_y(rad),0])
    r_mat = rotation_matrix(np.array([0,0,1]),pi/3)
    p2 = dot(r_mat,p1)
    p3 = dot(r_mat,p2)
    p4 = dot(r_mat,p3)
    p5 = dot(r_mat,p4)
    p6 = dot(r_mat,p5)
    return dwg.polygon(points = [(p[0],p[1]) for p in [p1,p2,p3,p4,p5,p6]],
        fill='none',
        stroke='black',
        stroke_width = 1)

def make_row(row_num):
    items = list()
    row_mod = (row_num%2)*(inside_rad+wall_thick)/2
    #print row_mod
    left_line = dwg.line(start = (-1*inside_rad-wall_thick,hex_y(inside_rad+wall_thick)),
        end = (-1*inside_rad-wall_thick,hex_y(inside_rad+wall_thick)*-1),
        fill='none',
        stroke='black',
        stroke_width = 1)
    tx = (row_mod+(inside_rad+wall_thick))*2
    ty = (inside_rad+wall_thick)*2*cos(pi/6)*row_num
    left_line.translate(tx=tx,ty=ty)
    items.append(left_line)
    cols_to_make = num_cols
    if not(row_num%2):
        cols_to_make += 1
    for hex_num in range(cols_to_make+1)[1:]:
        hex_group = dwg.g()
        hex = make_hex(inside_rad)
        tx = row_mod+(inside_rad+wall_thick)*hex_num
        ty = (inside_rad+wall_thick)*2*cos(pi/6)*row_num
        hex_group.add(hex)
        for h in make_hexes(inside_rad+shift,dowel_rad):
            hex_group.add(h)
        tx = (row_mod+(inside_rad+wall_thick)*hex_num)*2
        ty = (inside_rad+wall_thick)*2*cos(pi/6)*row_num


        right_line = dwg.line(start = (inside_rad+wall_thick,hex_y(inside_rad+wall_thick)),
            end = (inside_rad+wall_thick,hex_y(inside_rad+wall_thick)*-1),
            fill='none',
            stroke='black',
            stroke_width = 1)
        right_line.translate(tx=tx,ty=ty)
        hex_group.translate(tx=tx,ty=ty)
        items.append(hex_group)
        items.append(right_line)


    top_cut = tesselate_segment(inside_rad+wall_thick,cols_to_make,
        add_start = not((row_num+1)%2)-(row_num==1),
        add_end = 0)

    top_cut = dwg.polyline(top_cut,
        fill='none',
        stroke='black',
        stroke_width = 1)
    tx = (row_mod+(inside_rad+wall_thick))*2
    ty = (inside_rad+wall_thick)*2*cos(pi/6)*row_num
    top_cut.translate(tx = tx, ty = ty)
    items.append(top_cut)
    print row_num
    if row_num == num_rows:
        print 'here'
        row_mod = ((row_num+1)%2)*(inside_rad+wall_thick)/2
        bottom_cut = tesselate_segment(inside_rad+wall_thick,cols_to_make-1,
            add_start = 1,
            add_end = 1)
        bottom_cut = dwg.polyline(bottom_cut,
            fill='none',
            stroke='black',
            stroke_width = 1)
        tx = (inside_rad+wall_thick)*3
        ty = (inside_rad+wall_thick)*2*cos(pi/6)*(row_num+1)
        bottom_cut.translate(tx = tx, ty = ty)
        items.append(bottom_cut)
    return items

def tesselate_segment(rad,num, add_start= True ,add_end = True):
    p1 = np.array([rad,rad*tan(pi/6),0])
    r_mat = rotation_matrix(np.array([0,0,1]),pi/3)
    p2 = dot(r_mat,p1)
    p3 = dot(r_mat,p2)
    p4 = dot(r_mat,p3)
    first_seg = [(p[0],p[1]) for p in [p3,p2]]

    if add_start:
        pts = [(p3[0]-2*rad,p3[1])]
    else:
        pts = []
    pts.extend([(p4[0],p4[1])])
    for rep in range(num):
        pts.extend([(pt[0] + 2*(rep * rad),pt[1])
        for pt in first_seg])
    if add_end:
        pts.append((p3[0] + 2*(num * rad), p3[1]))
    return pts

def rotation_matrix(axis,theta):
    axis = axis/np.sqrt(np.dot(axis,axis))
    a = np.cos(theta/2)
    b,c,d = -axis*np.sin(theta/2)
    return np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
        [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
        [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]])

def make_cutout():
    pth = dwg.path(('m',(0,0)))
    seg_len = cutout_depth/(num_nubins)
    #print seg_len
    cut_pth = dwg.path(('m',(0,0)),
        fill='none',
        stroke='black',
        stroke_width = 1)
    for seg_num in range(num_nubins):
        cut_pth.push_arc((seg_len/2.0,0),0,(seg_len/4,nubin_d),angle_dir='-')
        cut_pth.push('l',(seg_len/2.0,0))
    cut_pth.push('l',(0,material_thick))
    for seg_num in range(num_nubins):
        cut_pth.push('l',(seg_len/2.0*-1,0))
        cut_pth.push_arc((seg_len/2.0*-1,0),0,(seg_len/4,nubin_d), angle_dir='-')
    return cut_pth


def add_cutouts(coord,orient):
    tvec = to_hex_coord(*coord)
    cutouts = [make_cutout() for x in range(2)]
    for i,cutout in enumerate(cutouts):
        cutout.translate(tx = tvec[0] + inside_rad, ty = tvec[1] - 0.5*material_thick)
        cutout.rotate(degrees((i*2 + orient)*pi/3),center = (inside_rad*-1,0.5*material_thick))
        dwg.add(cutout)


def make_rail():
    cutouts = [make_cutout() for x in range(4)]
    [c.translate(tx = 0,ty = sep*i) for i,c in enumerate(cutouts)]
    [c.rotate(180) for c in cutouts[::2]]
    [c.translate(tx = -1*wall_thick,ty = material_thick*-1) for c in cutouts[::2]]
    rail = dwg.g()
    outer_box = dwg.rect(size=(wall_thick,sep*(4-1)+material_thick + 2*front_thick),
        ry = wall_thick/10,
        fill='none',
        stroke='black',
        stroke_width = 1)
    outer_box.translate(tx = 0,ty = -1*front_thick)

    rail.add(outer_box)
    [rail.add(c) for c in cutouts]
    return rail

if __name__ == '__main__':
    dwg = svgwrite.Drawing('test.svg',size=(200*mm,200*mm),viewBox = '0 0 200 200')
    for row in range(num_rows+1)[1:]:
        items = make_row(row)
        for item in items:
            dwg.add(item)
    circ = dwg.circle(fill='none',
        stroke='black',
        stroke_width = 1,
        r=20)
    add_cutouts((1,1),0)
    add_cutouts((1,2),1)
    add_cutouts((1,3),2)
    add_cutouts((1,4),3)
    add_cutouts((2,1),4)
    add_cutouts((2,2),5)
    add_cutouts((2,3),1)
    add_cutouts((2,4),2)
    add_cutouts((2,5),3)
    dwg.add(make_rail())
    dwg.save()
