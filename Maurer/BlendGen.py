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
    if IsTLine(s): return SkipCharAny(s)
    if IsTElt(s):  return SkipCharAny(s)
    tf, val = IsNumber(s)
    if tf: return SkipNumber(s)

    if IsEmpty(s): return s
    assert 0

def GetTokens(s):
    w = []
    while not IsEmpty(s):
        w.append(Token(s))
        s = SkipToken(s)
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
    
    if len(tok):
        assert 0 # Not all tokens consumed

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

print('\n')
sss = '; ,1, 1, 1 ; ,0 ,0 ,0'
PrintTokens(sss)
print(GetRows(sss))

