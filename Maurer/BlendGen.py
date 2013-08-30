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
sss = """
; ,0 ,0 ,1 
; ,0 ,0 ,1 
; ,1 ,1 ,2 

; ,0 ,0 ,1 
; ,0 ,0 ,1 
; ,1 ,1 ,2 

; ,1 ,1 ,2 
; ,1 ,1 ,2 
; ,2 ,2 ,3"""
PrintTokens(sss)
print(GetRows(sss))

print(vLtoArray(GetRows(sss)))

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
def AddCube(origin):
    bpy.ops.mesh.primitive_cube_add(location=origin)
    return bpy.context.active_object

def AddCubeCol(origin, col):
    ob = AddCube(origin)
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

SetViewportShadeTextured()

ob = AddCubeCol([0,0,0], [0.3, 0.3, 0.3])
