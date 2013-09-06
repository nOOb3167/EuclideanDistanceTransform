sss0 = """
; ,0 ,0 ,1 
; ,0 ,0 ,1 
; ,1 ,1 ,2 

; ,0 ,0 ,1 
; ,0 ,0 ,1 
; ,1 ,1 ,2 

; ,1 ,1 ,2 
; ,1 ,1 ,2 
; ,2 ,2 ,3 
"""

sss1 = """
; ,0 ,0 ,1 ,4 ,9 ,16 ,25 ,36 ,49 ,64 
; ,0 ,0 ,1 ,4 ,9 ,16 ,25 ,36 ,49 ,57 
; ,1 ,1 ,2 ,5 ,10 ,17 ,26 ,37 ,43 ,50 
; ,4 ,4 ,5 ,8 ,13 ,20 ,29 ,33 ,38 ,45 
; ,9 ,9 ,10 ,13 ,18 ,25 ,27 ,30 ,35 ,42 
; ,16 ,16 ,17 ,20 ,25 ,25 ,26 ,29 ,34 ,41 
; ,25 ,25 ,26 ,29 ,27 ,26 ,27 ,30 ,35 ,42 
; ,36 ,36 ,37 ,33 ,30 ,29 ,30 ,33 ,38 ,45 
; ,49 ,49 ,43 ,38 ,35 ,34 ,35 ,38 ,43 ,50 
; ,64 ,57 ,50 ,45 ,42 ,41 ,42 ,45 ,50 ,57 

; ,0 ,0 ,1 ,4 ,9 ,16 ,25 ,36 ,49 ,57 
; ,0 ,0 ,1 ,4 ,9 ,16 ,25 ,36 ,41 ,48 
; ,1 ,1 ,2 ,5 ,10 ,17 ,26 ,29 ,34 ,41 
; ,4 ,4 ,5 ,8 ,13 ,20 ,21 ,24 ,29 ,36 
; ,9 ,9 ,10 ,13 ,18 ,17 ,18 ,21 ,26 ,33 
; ,16 ,16 ,17 ,20 ,17 ,16 ,17 ,20 ,25 ,32 
; ,25 ,25 ,26 ,21 ,18 ,17 ,18 ,21 ,26 ,33 
; ,36 ,36 ,29 ,24 ,21 ,20 ,21 ,24 ,29 ,36 
; ,49 ,41 ,34 ,29 ,26 ,25 ,26 ,29 ,34 ,41 
; ,57 ,48 ,41 ,36 ,33 ,32 ,33 ,36 ,41 ,48 

; ,1 ,1 ,2 ,5 ,10 ,17 ,26 ,37 ,43 ,50 
; ,1 ,1 ,2 ,5 ,10 ,17 ,26 ,29 ,34 ,41 
; ,2 ,2 ,3 ,6 ,11 ,18 ,19 ,22 ,27 ,34 
; ,5 ,5 ,6 ,9 ,14 ,13 ,14 ,17 ,22 ,29 
; ,10 ,10 ,11 ,14 ,11 ,10 ,11 ,14 ,19 ,26 
; ,17 ,17 ,18 ,13 ,10 ,9 ,10 ,13 ,18 ,25 
; ,26 ,26 ,19 ,14 ,11 ,10 ,11 ,14 ,19 ,26 
; ,37 ,29 ,22 ,17 ,14 ,13 ,14 ,17 ,22 ,29 
; ,43 ,34 ,27 ,22 ,19 ,18 ,19 ,22 ,27 ,34 
; ,50 ,41 ,34 ,29 ,26 ,25 ,26 ,29 ,34 ,41 

; ,4 ,4 ,5 ,8 ,13 ,20 ,29 ,33 ,38 ,45 
; ,4 ,4 ,5 ,8 ,13 ,20 ,21 ,24 ,29 ,36 
; ,5 ,5 ,6 ,9 ,14 ,13 ,14 ,17 ,22 ,29 
; ,8 ,8 ,9 ,12 ,9 ,8 ,9 ,12 ,17 ,24 
; ,13 ,13 ,14 ,9 ,6 ,5 ,6 ,9 ,14 ,21 
; ,20 ,20 ,13 ,8 ,5 ,4 ,5 ,8 ,13 ,20 
; ,29 ,21 ,14 ,9 ,6 ,5 ,6 ,9 ,14 ,21 
; ,33 ,24 ,17 ,12 ,9 ,8 ,9 ,12 ,17 ,24 
; ,38 ,29 ,22 ,17 ,14 ,13 ,14 ,17 ,22 ,29 
; ,45 ,36 ,29 ,24 ,21 ,20 ,21 ,24 ,29 ,36 

; ,9 ,9 ,10 ,13 ,18 ,25 ,27 ,30 ,35 ,42 
; ,9 ,9 ,10 ,13 ,18 ,17 ,18 ,21 ,26 ,33 
; ,10 ,10 ,11 ,14 ,11 ,10 ,11 ,14 ,19 ,26 
; ,13 ,13 ,14 ,9 ,6 ,5 ,6 ,9 ,14 ,21 
; ,18 ,18 ,11 ,6 ,3 ,2 ,3 ,6 ,11 ,18 
; ,25 ,17 ,10 ,5 ,2 ,1 ,2 ,5 ,10 ,17 
; ,27 ,18 ,11 ,6 ,3 ,2 ,3 ,6 ,11 ,18 
; ,30 ,21 ,14 ,9 ,6 ,5 ,6 ,9 ,14 ,21 
; ,35 ,26 ,19 ,14 ,11 ,10 ,11 ,14 ,19 ,26 
; ,42 ,33 ,26 ,21 ,18 ,17 ,18 ,21 ,26 ,33 

; ,9 ,9 ,10 ,13 ,18 ,25 ,26 ,29 ,34 ,41 
; ,9 ,9 ,10 ,13 ,17 ,16 ,17 ,20 ,25 ,32 
; ,10 ,10 ,11 ,13 ,10 ,9 ,10 ,13 ,18 ,25 
; ,13 ,13 ,13 ,8 ,5 ,4 ,5 ,8 ,13 ,20 
; ,18 ,17 ,10 ,5 ,2 ,1 ,2 ,5 ,10 ,17 
; ,25 ,16 ,9 ,4 ,1 ,0 ,1 ,4 ,9 ,16 
; ,26 ,17 ,10 ,5 ,2 ,1 ,2 ,5 ,10 ,17 
; ,29 ,20 ,13 ,8 ,5 ,4 ,5 ,8 ,13 ,20 
; ,34 ,25 ,18 ,13 ,10 ,9 ,10 ,13 ,18 ,25 
; ,41 ,32 ,25 ,20 ,17 ,16 ,17 ,20 ,25 ,32 

; ,4 ,4 ,5 ,8 ,13 ,20 ,27 ,30 ,35 ,42 
; ,4 ,4 ,5 ,8 ,13 ,17 ,18 ,21 ,26 ,33 
; ,5 ,5 ,6 ,9 ,11 ,10 ,11 ,14 ,19 ,26 
; ,8 ,8 ,9 ,9 ,6 ,5 ,6 ,9 ,14 ,21 
; ,13 ,13 ,11 ,6 ,3 ,2 ,3 ,6 ,11 ,18 
; ,20 ,17 ,10 ,5 ,2 ,1 ,2 ,5 ,10 ,17 
; ,27 ,18 ,11 ,6 ,3 ,2 ,3 ,6 ,11 ,18 
; ,30 ,21 ,14 ,9 ,6 ,5 ,6 ,9 ,14 ,21 
; ,35 ,26 ,19 ,14 ,11 ,10 ,11 ,14 ,19 ,26 
; ,42 ,33 ,26 ,21 ,18 ,17 ,18 ,21 ,26 ,33 

; ,1 ,1 ,2 ,5 ,10 ,17 ,26 ,33 ,38 ,45 
; ,1 ,1 ,2 ,5 ,10 ,17 ,21 ,24 ,29 ,36 
; ,2 ,2 ,3 ,6 ,11 ,13 ,14 ,17 ,22 ,29 
; ,5 ,5 ,6 ,9 ,9 ,8 ,9 ,12 ,17 ,24 
; ,10 ,10 ,11 ,9 ,6 ,5 ,6 ,9 ,14 ,21 
; ,17 ,17 ,13 ,8 ,5 ,4 ,5 ,8 ,13 ,20 
; ,26 ,21 ,14 ,9 ,6 ,5 ,6 ,9 ,14 ,21 
; ,33 ,24 ,17 ,12 ,9 ,8 ,9 ,12 ,17 ,24 
; ,38 ,29 ,22 ,17 ,14 ,13 ,14 ,17 ,22 ,29 
; ,45 ,36 ,29 ,24 ,21 ,20 ,21 ,24 ,29 ,36 

; ,0 ,0 ,1 ,4 ,9 ,16 ,25 ,36 ,43 ,50 
; ,0 ,0 ,1 ,4 ,9 ,16 ,25 ,29 ,34 ,41 
; ,1 ,1 ,2 ,5 ,10 ,17 ,19 ,22 ,27 ,34 
; ,4 ,4 ,5 ,8 ,13 ,13 ,14 ,17 ,22 ,29 
; ,9 ,9 ,10 ,13 ,11 ,10 ,11 ,14 ,19 ,26 
; ,16 ,16 ,17 ,13 ,10 ,9 ,10 ,13 ,18 ,25 
; ,25 ,25 ,19 ,14 ,11 ,10 ,11 ,14 ,19 ,26 
; ,36 ,29 ,22 ,17 ,14 ,13 ,14 ,17 ,22 ,29 
; ,43 ,34 ,27 ,22 ,19 ,18 ,19 ,22 ,27 ,34 
; ,50 ,41 ,34 ,29 ,26 ,25 ,26 ,29 ,34 ,41 

; ,0 ,0 ,1 ,4 ,9 ,16 ,25 ,36 ,49 ,57 
; ,0 ,0 ,1 ,4 ,9 ,16 ,25 ,36 ,41 ,48 
; ,1 ,1 ,2 ,5 ,10 ,17 ,26 ,29 ,34 ,41 
; ,4 ,4 ,5 ,8 ,13 ,20 ,21 ,24 ,29 ,36 
; ,9 ,9 ,10 ,13 ,18 ,17 ,18 ,21 ,26 ,33 
; ,16 ,16 ,17 ,20 ,17 ,16 ,17 ,20 ,25 ,32 
; ,25 ,25 ,26 ,21 ,18 ,17 ,18 ,21 ,26 ,33 
; ,36 ,36 ,29 ,24 ,21 ,20 ,21 ,24 ,29 ,36 
; ,49 ,41 ,34 ,29 ,26 ,25 ,26 ,29 ,34 ,41 
; ,57 ,48 ,41 ,36 ,33 ,32 ,33 ,36 ,41 ,48 
"""

class E:
    TLINE = 0
    TELT  = 1
    TNUM  = 2
    TEOF  = 3

def sw(s):
    while len(s) and IsCharWs(s[0]):
        s = s[1:]
    return s

def IsChar(s, c):
    return len(s) and s[0] == c

def IsCharWs(s):
    return IsChar(s, ' ') or IsChar(s, '\t') or IsChar(s, '\n')

def IsCharNum(s):
    return any([IsChar(s, i) for i in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']])

def IsNumber(s):
    t = s
    while IsCharNum(t):
        t = t[1:]
    
    if len(t) != 0: sliceXLen = -len(t)
    else:           sliceXLen = None
    
    try:    return [True, int(s[:sliceXLen])]
    except: return [False, None]

def IsEmpty(s):
    return not len(s)

def SkipCharAny(s):
    if not len(s): return s
    return s[1:]

def SkipNumber(s):
    while IsCharNum(s):
        s = s[1:]
    return s

def IsTLine(s):
    return IsChar(s, ';')

def IsTElt(s):
    return IsChar(s, ',')

def Token(s):
    s = sw(s)
    if IsTLine(s): return [E.TLINE, ';']
    if IsTElt(s):  return [E.TELT, ',']
    tf, val = IsNumber(s)
    if tf: return [E.TNUM, val]

    if IsEmpty(s): return [E.TEOF, None]
    assert 0

def SkipToken(s):
    s = sw(s)
    
    if IsEmpty(s): return s
    
    if IsTLine(s): return SkipCharAny(s)
    if IsTElt(s):  return SkipCharAny(s)

    tf, val = IsNumber(s)
    if tf: return SkipNumber(s)

    assert 0

def GetTokens(s):
    w = []
    while True:
        t, ex = Token(s)
        w.append(Token(s))
        s = SkipToken(s)
        
        if t == E.TEOF:
            break
    return w

def PrintTokens(s):
    for t, ex in GetTokens(s):
        print('Token', t, 'Extra', ex)
    print('End')

def ParseElts(vL, tok):
    w = []
    
    while len(tok):
        t, ex = tok[0]

        if t != E.TELT:
            break
        else:
            tok = tok[1:]
            # A ',' requires a following 'NUM'
            assert len(tok)
            t, ex = tok[0]
            assert t == E.TNUM
            w.append(ex)
            tok = tok[1:]
    
    vL.append(w)
    return tok

def ParseRows(tok):
    vL = []

    while len(tok):
        t, ex = tok[0]
        
        if t != E.TLINE:
            break
        else:
            tok = tok[1:]
            tok = ParseElts(vL, tok)
    
    # Only EOF should be remaining. If > 1 not all tokens were consumed.
    assert len(tok) == 1
    assert tok[0][0] == E.TEOF and tok[0][1] == None

    return vL

def CheckSizesRowmatch(vL):
    assert len(vL)
    for i in vL:
        assert len(vL[0]) == len(i)

def GetRows(s):
    tok = GetTokens(s)
    vL = ParseRows(tok)
    CheckSizesRowmatch(vL)
    return vL

def GetUniformDims3(vL):
    assert len(vL)
    colLen = len(vL[0])
    rowsNeeded = colLen * colLen
    
    assert len(vL) == rowsNeeded
    
    return [colLen, colLen, colLen]

def vLtoArray(vL):
    dims = GetUniformDims3(vL)
    assert dims[0] == dims[1] and dims[0] == dims[2]
    f = MakeUniform3D(0, dims[0])
    
    for s in range(dims[2]):
        for r in range(dims[1]):
            for c in range(dims[0]):
                f[s][r][c] = vL[(dims[1] * s) + r][c]
    
    return f

######
def Make3D(init, z, y, x):
        r = []
        for k in range(z):
                r.append([])
                for j in range(y):
                        r[k].append([])
                        for i in range(x):
                                r[k][j].append(init)
        return r

def MakeUniform3D(init, n):
        return Make3D(init, n, n, n)
######

print('\n')
#sss = '; ,1, 1, 1 ; ,0 ,0 ,0 ; ,0 ,0 ,0'
#PrintTokens(sss0)
#print(GetRows(sss0))
#print(vLtoArray(GetRows(sss0)))

import math
import bpy

def makeMaterial(name, diffuse, specular, alpha):
    mat = bpy.data.materials.new(name)
    mat.diffuse_color = diffuse
    mat.diffuse_shader = 'LAMBERT' 
    mat.diffuse_intensity = 1.0 
    mat.specular_color = specular
    mat.specular_shader = 'COOKTORR'
    mat.specular_intensity = 0.5
    mat.alpha = alpha
    mat.ambient = 1
    return mat

def setMaterial(ob, mat):
    me = ob.data
    me.materials.append(mat)

# r = Object, r.data = datablock
def AddCubeSiz(origin, siz):
    bpy.ops.object.select_all(action='DESELECT')
    bpy.ops.mesh.primitive_cube_add(location=origin)
    bpy.ops.transform.resize(value=siz)
    return bpy.context.active_object

def AddCube(origin):
    return AddCubeSiz(origin, [0.5, 0.5, 0.5])

def AddCubeSizCol(origin, siz, col):
    ob = AddCubeSiz(origin, siz)
    me = ob.data
    me.vertex_colors.new(name='Col')
    cmap = me.vertex_colors['Col']
    for i in cmap.data:
        i.color = col        
    return ob

def SetViewportShadeTextured():
    for scrn in bpy.data.screens:
        if scrn.name == 'Default' or scrn.name == 'Scripting':
            for area in scrn.areas:
                if area.type == 'VIEW_3D':
                    for space in area.spaces:
                        if space.type == 'VIEW_3D':
                            space.viewport_shade = 'TEXTURED'

def NukeMeshObjs():
    for obname in [item.name for item in bpy.data.objects if item.type == "MESH"]:
        bpy.data.objects[obname].select = True
        
    bpy.ops.object.delete()
    
    for item in bpy.data.meshes:
        bpy.data.meshes.remove(item)

def GetRDims(f):
    return [len(f), len(f[0]), len(f[0][0])]

def GetRColorScalingFactor(f):
    dims = GetRDims(f)
    distMaxMax = dims[0] * dims[0] + dims[1] * dims[1] + dims[2] * dims[2]
    return 1.0 / math.sqrt(distMaxMax)

def MakeObjGrid(f):
    dims = GetRDims(f)
    csf = GetRColorScalingFactor(f)
    for k in range(dims[2]):
        for j in range(dims[1]):
            for i in range(dims[0]):
                v = f[i][j][k]
                c = [v*csf, v*csf, v*csf]
                AddCubeSizCol([i, j, k], [0.1, 0.5, 0.5], c)

NukeMeshObjs()

SetViewportShadeTextured()

#ob = AddCubeCol([0,0,0], [0.3, 0.3, 0.3])

aRows = vLtoArray(GetRows(sss1))
MakeObjGrid(aRows)
