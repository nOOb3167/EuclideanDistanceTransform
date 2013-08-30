class D:
    TLINE = 0
    TELT  = 1
    
    def __init__(self, s):
        self.s = s
        self.p = 0

    def Eof(self):
        return self.p >= len(self.s)
        
    def OptSkipWs(self):
        while not self.Eof() and IsCharWs(self.PeekChar()):
            self.p += 1
            
    def PeekChar(self):
        assert not self.Eof()
        return self.s[self.p]
    
    def GetChar(self):
        assert not self.Eof()
        c = self.s[self.p]
        self.p += 1
        return c
    
    def CToken(self):
        gc = self.GetChar()
        print('gcx',gc)
        if gc == ';': return self.TLINE
        if gc == ',': return self.TELT
        assert False

    def Unget(self):
        assert self.p > 0
        self.p -= 1
        
    def ExNum(self):
        assert not self.Eof()
        
        tmpP = self.p
        while not tmpP >= len(self.s) and IsCharNum(self.s[tmpP]):
            tmpP += 1
            
        w = int(self.s[self.p : tmpP])
        self.p = tmpP
        return w
        
    def ExElt(self):
        self.OptSkipWs()
        w = self.ExNum()
        self.OptSkipWs()
        return w

class E:
    TLINE = 0
    TELT  = 1
    TNUM  = 2
    TEOF  = 3

    def __init__(self, s):
        self.s = s
        self.p = 0

    def Token(self):
        pass

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

def PrintTokens(s):
    while not IsEmpty(s):
        t, ex = Token(s)
        s = SkipToken(s)
        print('Token', t, 'Extra', ex)
    print('End')

print('\n')
PrintTokens('; ,1, 1, 1 ; ,0 ,0 ,0')

#print(IsNumber('   \n\n\t 123 a bcdef'))
#print(SkipNumber('   \n\n\t 123 a bcdef'))

def xstate(d):
    print('|', d.s[:d.p], '| @ |', d.s[d.p:], '|')

def IsCharWs(c):
    return c == ' ' or c == '\t' or c == '\n'

def IsCharNum(c):
    return c in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']

def CheckSizes1(vvLines):
#    print(vvLines)
#    for i in range(len(vvLines)):
#            assert len(vvLines[i]) == 1
    pass

def ParseRow(d):
    vLine = []
    
    d.OptSkipWs()
    
    while not d.Eof() and d.CToken() == D.TELT:
        xstate(d)
        arityoneElt = d.ExElt()
        vLine.append(arityoneElt)
        
    if not d.Eof():
        d.Unget()
        
    return vLine

def ParseRows(d):
    vvLines = []
    
    d.OptSkipWs()
    
    while not d.Eof() and d.CToken() == D.TLINE:
        vvLines.append(ParseRow(d))
    
    if not d.Eof():
        d.Unget()
        
    return vvLines
    
def GetRows(s):
    d = D(s)
    vvLines = ParseRows(d)
    assert d.Eof()
    CheckSizes1(vvLines)
    return vvLines

#v = GetRows('; ,1, 1, 1 ; 0, 0, 0')
#print(v)
