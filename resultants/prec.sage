load "/home/xavier/git/padicprec/sage/default.sage"

def gcd(A,B):
    degree = A.degree()
    #print poly_to_string(A, degree, dec=1)
    precs = [ ] # (A.degree()+1, precision(A.monic())) ]
    while B != 0:
        precs.append( (B.degree(), precision(B.monic())) )
        #print poly_to_string(B, degree, dec=1)
        A, B = B, rem(A,B)
    return A.monic(), precs

def subres(A,B):
    degree = A.degree()
    #print poly_to_string(A, degree, dec=1)
    precs = [ ] # (A.degree()+1, precision(A.monic())) ]
    while B != 0:
        precs.append( (B.degree(), precision(B)) )
        #print poly_to_string(B, degree, dec=1)
        A, B = B, rem(B.leading_coefficient()^2*A,B) / A.leading_coefficient()^2
    return A.monic(), precs

def padic_to_string(x, dec=0):
    p = x.parent().prime()
    prec = x.precision_absolute()
    if prec < 0: raise NotImplementedError
    val = min(0, x.valuation())
    x = ZZ(x >> val); s = '';
    for i in range(prec-val):
        s = str(x % p) + s
        if i == -val-1:
            s = ",\\!" + s
        x = x // 2
    s = "\\ldots %s" % s
    if dec > 0:
        phantom = '';
        if val == 0: phantom = ",\\!"
        if -val < dec: phantom += (dec+val) * "0"
        if phantom != "": s += "\hphantom{%s}" % phantom
    return s

def poly_to_string(P, until=None, dec=0):
    deg = P.degree()
    val = P.leading_coefficient().valuation()
    if until is None: until = deg
    ans = "p^{%s} & " % val;
    for i in range(until,-1,-1):
        if i > deg: s = ''
        else: s = padic_to_string(P[i] >> val, dec=dec)
        ans += s + ' & '
    ans = ans[:-2] + "\\medskip \\\\";
    return ans
    


def tikz_output(points, filename="output.tex"):
    absc = [ x for (x,_,_) in points ]
    ordo = [ y for (_,y,_) in points ]
    xmin = 0; xmax = max(absc); dx = xmax - xmin
    ymin = 0; ymax = max(ordo); dy = ymax - ymin

    xscale = RR(9/dx)
    yscale = RR(6/dy)
    if dx < 25:
      xellipse = 0.05 / xscale
      yellipse = 0.05 / yscale
    elif dx < 50:
      xellipse = 0.02 / xscale
      yellipse = 0.02 / yscale
    else:
      xellipse = 0.01 / xscale
      yellipse = 0.01 / yscale

    handle = open(filename, "w")
    handle.write("\\begin{tikzpicture}[xscale=%s,yscale=%s]\n" % (xscale, yscale))
    handle.write("\\draw[->] (%s,%s)--(%s,%s);\n" % (xmin-dx/15,ymin-dy/30,xmin-dx/15,ymax+dy/30))
    handle.write("\\draw[->] (%s,%s)--(%s,%s);\n" % (xmin-dx/30,ymin-dy/15,xmax+dx/30,ymin-dy/15))

    y0 = ymin - dy/15 - 1/(20*yscale)
    y1 = ymin - dy/15 + 1/(20*yscale)
    if dx < 10: step = 1
    elif dx < 100: step = 10
    elif dx < 200: step = 20
    elif dx < 500: step = 50
    else: step = 100
    for x in range(xmin,xmax+1):
        if x % step != 0: continue
        handle.write("\\draw (%s,%s)--(%s,%s);\n" % (x,y0,x,y1))
        handle.write("\\node[scale=0.5,below] at (%s,%s) { $%s$ };\n" % (x,y0,x))
    x0 = xmin - dx/15 - 1/(20*xscale)
    x1 = xmin - dx/15 + 1/(20*xscale)
    if dy < 10: step = 1
    elif dy < 100: step = 10
    elif dy < 200: step = 20
    elif dy < 500: step = 50
    else: step = 100
    for y in range(ymin,ymax+1):
        if y % step != 0: continue
        handle.write("\\draw (%s,%s)--(%s,%s);\n" % (x0,y,x1,y))
        handle.write("\\node[scale=0.5,left] at (%s,%s) { $%s$ };\n" % (x0,y,y))

    for point in points:
        handle.write("\\draw[%s,fill] (%s,%s) ellipse (%s and %s);\n" % (point[2], point[0], point[1], xellipse, yellipse))
        first = False
    handle.write("\\end{tikzpicture}\n")
    handle.close()


def blah(p=2, degree=20, prec=50, A=None, B=None, filename="output.tex"):
  R = Zp(p, prec=prec)
  S.<x> = PolynomialRing(R)
  if A is None:
      A = random_monic(S, degree)
  else:
      A = S(A)
  if B is None:
      B = random_monic(S, degree)
  else:
      B = S(B)

  k = R.residue_field()
  Sk.<x> = PolynomialRing(k)
  Abar = Sk(A); Bbar = Sk(B)
  if Abar.resultant(Bbar) == 0:
      return

  D, precs = subres(A,B)
  #tikz_output(precs, filename=filename)
  return precs, A, B


def blah2(val, shift=15, ymax=45, filename="output.tex"):
    degree = len(val)

    xmin = 0; xmax = degree+1; dx = xmax - xmin
    ymin = 0; dy = ymax - ymin

    xscale = RR(9/dx)
    yscale = RR(6/dy)
    if dx < 25:
      xellipse = 0.05 / xscale
      yellipse = 0.05 / yscale
    elif dx < 50:
      xellipse = 0.02 / xscale
      yellipse = 0.02 / yscale
    else:
      xellipse = 0.01 / xscale
      yellipse = 0.01 / yscale

    handle = open(filename, "w")
    handle.write("\\begin{tikzpicture}[xscale=%s,yscale=%s]\n" % (xscale, yscale))
    handle.write("\\draw[->] (%s,%s)--(%s,%s);\n" % (xmin-dx/15,ymin-dy/30,xmin-dx/15,ymax+dy/30))
    handle.write("\\draw[->] (%s,%s)--(%s,%s);\n" % (xmin-dx/30,ymin-dy/15,xmax+dx/30,ymin-dy/15))

    y0 = ymin - dy/15 - 1/(20*yscale)
    y1 = ymin - dy/15 + 1/(20*yscale)
    if dx < 10: step = 1
    elif dx < 100: step = 10
    elif dx < 200: step = 20
    elif dx < 500: step = 50
    else: step = 100
    for x in range(xmin,xmax+1):
        if x % step != 0: continue
        handle.write("\\draw (%s,%s)--(%s,%s);\n" % (x,y0,x,y1))
        handle.write("\\node[scale=0.5,below] at (%s,%s) { $%s$ };\n" % (x,y0,x))
    x0 = xmin - dx/15 - 1/(20*xscale)
    x1 = xmin - dx/15 + 1/(20*xscale)
    if dy < 10: step = 1
    elif dy < 100: step = 10
    elif dy < 200: step = 20
    elif dy < 500: step = 50
    else: step = 100
    for y in range(ymin,ymax+1):
        if y % step != 0: continue
        handle.write("\\draw (%s,%s)--(%s,%s);\n" % (x0,y,x1,y))
        handle.write("\\node[scale=0.5,left] at (%s,%s) { $%s$ };\n" % (x0,y,y))

    y_violet = [ ]
    for _ in range(degree+2):
        y_violet.append([])

    val = [ 0, 0 ] + val; val.reverse()
    handle.write("\\draw[rouge,fill] (%s,%s) ellipse (%s and %s);\n" % (degree+1, shift, xellipse, yellipse))
    frame = 1
    for j in range(degree, 0, -1):
        frame_end = frame + 1
        if val[j] > 0: frame_end += 2
        if val[j+1] > 0: frame_end += 1
        handle.write("\\only<%s-%s>{\n" % (frame, frame_end))
        handle.write("\\node[left] at (%s,%s) { \\ph $j =$ };\n" % (1.2/xscale, ymax - 0.5/yscale))
        handle.write("\\node[right] at (%s,%s) { \\ph $%s$ };\n" % (1.05/xscale, ymax - 0.5/yscale, j))
        handle.write("\\node[left] at (%s,%s) { \\ph $V_j =$ };\n" % (1.2/xscale, ymax - 1/yscale))
        handle.write("\\node[right] at (%s,%s) { \\ph $%s$ };\n" % (1.05/xscale, ymax - 1/yscale, val[j]))
        handle.write("\\node[left] at (%s,%s) { \\ph $V_{j+1} =$ };\n" % (1.2/xscale, ymax - 1.5/yscale))
        handle.write("\\node[right] at (%s,%s) { \\ph $%s$ };\n" % (1.05/xscale, ymax - 1.5/yscale, val[j+1]))
        handle.write("}\n")
        if val[j] > 0:
            y_violet[j].append([ shift + 2*val[j+1], frame, frame+1 ])
            y_violet[j+1].append([ shift + 2*val[j+1], frame, frame+1 ])
            handle.write("\\only<%s-%s>{\n" % (frame, frame+1))
            handle.write("\\draw[violet,fill] (%s,%s) ellipse (%s and %s);\n" % (j+1, shift + 2*val[j+1], xellipse, yellipse))
            handle.write("\\draw[violet,fill] (%s,%s) ellipse (%s and %s);\n" % (j, shift + 2*val[j+1], xellipse, yellipse))
            handle.write("}\n")
            frame += 1
            y_violet[j].append([ shift + 2*val[j] + 2*val[j+1], frame, frame_end ])
            y_violet[j+1].append([ shift + 2*val[j] + 2*val[j+1], frame, frame_end ])
            handle.write("\\only<%s-%s>{\n" % (frame, frame_end))
            handle.write("\\draw[violet,fill] (%s,%s) ellipse (%s and %s);\n" % (j+1, shift + 2*val[j] + 2*val[j+1], xellipse, yellipse))
            handle.write("\\draw[violet,fill] (%s,%s) ellipse (%s and %s);\n" % (j, shift + 2*val[j] + 2*val[j+1], xellipse, yellipse))
            handle.write("}\n")
            handle.write("\\only<%s>{\n" % frame)
            handle.write("\\draw[violet,->] (%s,%s) -- (%s,%s);\n" % (j+1, shift + 2*val[j+1] + 0.5, j+1, shift + 2*val[j] + 2*val[j+1] - 0.5))
            handle.write("\\draw[violet,->] (%s,%s) -- (%s,%s);\n" % (j, shift + 2*val[j+1] + 0.5, j, shift + 2*val[j] + 2*val[j+1] - 0.5))
            handle.write("}\n")
            frame += 2
        else:
            y_violet[j].append([ shift + 2*val[j] + 2*val[j+1], frame, frame_end ])
            y_violet[j+1].append([ shift + 2*val[j] + 2*val[j+1], frame, frame_end ])
            handle.write("\\only<%s-%s>{\n" % (frame, frame_end))
            handle.write("\\draw[violet,fill] (%s,%s) ellipse (%s and %s);\n" % (j+1, shift + 2*val[j+1], xellipse, yellipse))
            handle.write("\\draw[violet,fill] (%s,%s) ellipse (%s and %s);\n" % (j, shift + 2*val[j+1], xellipse, yellipse))
            handle.write("}\n")
            frame += 1
        y_violet[j-1].append([ shift + 2*val[j], frame, frame_end ])
        handle.write("\\only<%s-%s>{\n" % (frame, frame_end))
        handle.write("\\draw[violet,fill] (%s,%s) ellipse (%s and %s);\n" % (j-1, shift + 2*val[j], xellipse, yellipse))
        handle.write("}\n")
        frame += 1
        if val[j+1] > 0:
            y_violet[j].append([ shift + 2*val[j] + 2*val[j+1], frame, frame_end ])
            handle.write("\\only<%s-%s>{\n" % (frame, frame_end))
            handle.write("\\draw[violet,fill] (%s,%s) ellipse (%s and %s);\n" % (j, shift + 2*val[j], xellipse, yellipse))
            handle.write("}\n")
            handle.write("\\only<%s>{\n" % frame)
            handle.write("\\draw[violet,->] (%s,%s) -- (%s,%s);\n" % (j, shift + 2*val[j] + 2*val[j+1] - 0.5, j, shift + 2*val[j] + 0.5))
            handle.write("}\n")
            frame += 1
        handle.write("\\only<%s->{\n" % frame)
        if j == degree:
            handle.write("\\draw[rouge,fill] (%s,%s) ellipse (%s and %s);\n" % (j, shift + 2*val[j], xellipse, yellipse))
        #else:
        #    handle.write("\\draw[vert,fill] (%s,%s) ellipse (%s and %s);\n" % (j, shift + 2*val[j], xellipse, yellipse))
        handle.write("}\n")

    for j in range(degree):
        start = None; height = 0
        for (y, f1, f2) in y_violet[j]:
            if height > y: continue
            if height > 0 and f1-1 > start:
                handle.write("\\only<%s-%s>{ " % (start, f1-1))
                handle.write("\\draw[violet!50,fill] (%s,%s) ellipse (%s and %s);" % (j, height, xellipse, yellipse))
                handle.write(" }\n")
            start = f2 + 1
            height = max(height,y)
        if height > 0:
            handle.write("\\only<%s-%s>{ " % (start, frame))
            handle.write("\\draw[violet!50,fill] (%s,%s) ellipse (%s and %s);" % (j, height, xellipse, yellipse))
            handle.write(" }\n")
    print y_violet[19]

    handle.write("\\end{tikzpicture}\n")
    handle.close()


def blah3(val, shift=15, ymax=45, filename="output.tex"):
    degree = len(val)
    pts = [ ]
    val = [ 0, 0 ] + val; val.reverse()
    for j in range(degree-1, -1, -1):
        if j > 0:
            height = shift + 2*(val[j] + max(val[j+1], val[j-1]))
        else:
            height = shift + 2*(val[j] + val[j+1])
        pts.append((j, height, 'vert2'))
    return pts
