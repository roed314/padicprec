from sage.misc.misc import walltime, cputime

from sage.rings.polynomial.polynomial_ring import PolynomialRing_general
from sage.rings.power_series_ring import PowerSeriesRing_generic
from sage.matrix.matrix_space import MatrixSpace
import io

import pygments
lexer = pygments.lexers.PythonLexer()
html_formatter = pygments.formatters.HtmlFormatter()
# latex_formatter = pygments.formatters.LatexFormatter()

from nbconvert.filters.markdown import markdown2html

MAXINT = 2^1000


randints = [ ]
def random_element(ring, prec=None, degree=None, integral=False):
    global random_seed
    random_seed += 1
    for _ in range(len(randints), random_seed+2):
        randints.append(ZZ.random_element(MAXINT))
    if isinstance(ring, PolynomialRing_general):
        if degree is None:
            degree = 2
        R = ring.base_ring()
        coeffs = [ random_element(R, prec=prec, integral=integral) for _ in range(degree+1) ]
        return ring(coeffs).change_ring(R.fraction_field())
    elif isinstance(ring, PowerSeriesRing_generic):
        if degree is None:
            degree = ring.default_prec()
        R = ring.base_ring()
        coeffs = [ random_element(R, prec=prec, integral=integral) for _ in range(degree+1) ]
        return ring(coeffs)
    elif isinstance(ring, MatrixSpace):
        size = ring.ncols() * ring.nrows()
        R = ring.base_ring()
        coeffs = [ random_element(R, prec=prec, degree=degree, integral=integral) for _ in range(size) ]
        return ring(coeffs)
    else:
        if not integral and ring.is_field():
            r = randints[random_seed] - MAXINT/2
            random_seed += 1
            val = floor(2*MAXINT/(5*r))
        else:
            val = 0
        if prec is None:
            return ring(randints[random_seed]) << val
        else:
            return (ring(randints[random_seed]).add_bigoh(prec) << val)


def my_exec(code, ring=Zp, field=Qp):
    global random_seed
    random_seed = -1;
    Zp = ring
    Qp = field
    Out = [ ]
    Tme = [ ]
    for cmd in code:
        Out.append(None)
        Tme.append(None)
        length = len(cmd)
        if length == 0: continue
        if cmd[0] == "%% markdown\n": continue
        cmd_exec = ""
        for n in range(length-1):
            cmd_exec += preparse(cmd[n])
        cmd_eval = preparse(cmd[-1])
        tme = cputime()
        exec(cmd_exec)
        try:
            Out[-1] = eval(cmd_eval)
        except SyntaxError:
            exec(cmd_eval)
        except Exception as e:
            Out[-1] = e.__class__.__name__ + ": " + str(e)
        Tme[-1] = cputime(tme)
    return Out, Tme
    

def readfile(filename):
    fh = io.open(filename)
    code = [ ]; cmd = [ ]
    while True:
        line = fh.readline()
        if not line: break
        if line == "\n":
            if len(cmd) > 0: code.append(cmd)
            cmd = [ ]
        else:
            if line[0] == " ":
                cmd[-1] += line
            else:
                cmd.append(line)
    if len(cmd) > 0: code.append(cmd)
    return code


def to_txt(filename, outfile=None):
    fh = None
    if outfile is not None:
        fh = io.open(outfile, "w")
    code = readfile(filename)
    OutCR = my_exec(code, ZpCR, QpCR)
    OutFP = my_exec(code, ZpFP, QpFP)
    OutLC = my_exec(code, ZpLC, QpLC)
    OutLF = my_exec(code, ZpLF, QpLF)
    s = ""
    for n in range(len(code)):
        for line in code[n]:
            lines = line.split("\n");
            s += "sage: " + lines[0] + "\n"
            for i in range(1,len(lines)-1):
                s += "....: " + lines[i] + "\n"
        if OutCR[n] is not None:
            s += "ZpCR:\n" + str(OutCR[n]) + "\n"
        if OutFP[n] is not None:
            s += "ZpFP:\n" + str(OutFP[n]) + "\n"
        if OutLC[n] is not None:
            s += "ZpLC:\n" + str(OutLC[n]) + "\n"
        if OutLF[n] is not None:
            s += "ZpLF:\n" + str(OutLF[n]) + "\n"
        s += "\n"

    if fh:
        fh.write(s)
        fh.close()
    return s


def htmlspecialchars(text):
    return (
        text.replace("&", "&amp;").
        replace('"', "&quot;").
        replace("<", "&lt;").
        replace(">", "&gt;")
    )
def to_html(filename, outfile=None):
    fh = None
    if outfile is not None:
        fh = io.open(outfile, "w")
    code = readfile(filename)
    OutCR, TmeCR = my_exec(code, ZpCR, QpCR)
    OutFP, TmeFP = my_exec(code, ZpFP, QpFP)
    OutLC, TmeLC = my_exec(code, ZpLC, QpLC)
    OutLF, TmeLF = my_exec(code, ZpLF, QpLF)

    fhead = io.open("head.html")
    if fhead:
        s = fhead.read()
        fhead.close()
    else:
        s = ""

    nn = 0
    for n in range(len(code)):
        s += " <div class=\"cell border-box-sizing code_cell rendered\">\n"

        # Input
        s += "  <div class=\"input\">\n"
        if code[n][0] == "%% markdown\n":
            no = ""
            s += "    <div class=\"prompt input_prompt\"></div>\n"
            s += "    <div class=\"inner_cell\">\n"
            s += "      <div class=\"text_cell_render border-box-sizing rendered_html\">\n"
            s += markdown2html("".join(code[n][1:])) + "\n"
        else:
            nn += 1
            no = "&nbsp;[" + str(nn) + "]:"
            s += "    <div class=\"prompt input_prompt\">In" + no + "</div>\n"
            s += "    <div class=\"inner_cell\">\n"
            s += "      <div class=\"input_area\">\n"
            s += pygments.highlight("".join(code[n]), lexer, html_formatter) + "\n"
        s += "      </div>\n";
        s += "    </div>\n";
        s += "  </div>\n";

        # Output
        if OutCR[n] is not None:
            s += "  <div class=\"output_wrapper\">\n"
            s += "    <div class=\"output\">\n"
            s += "      <div class=\"output_area\">\n"
            s += "        <div class=\"prompt output_prompt\">ZpCR" + no + "</div>\n"
            s += "        <div class=\"output_text output_subarea output_execute_result\">\n";
            s += "          <div class=\"timings\">%.3f&nbsp;s</div>\n" % TmeCR[n]
            s += "          <pre>" + htmlspecialchars(str(OutCR[n])) + "</pre>\n"
            s += "        </div>\n";
            s += "      </div>\n";
            s += "    </div>\n";
            s += "  </div>\n";
        if OutFP[n] is not None:
            s += "  <div class=\"output_wrapper\">\n"
            s += "    <div class=\"output\">\n"
            s += "      <div class=\"output_area\">\n"
            s += "        <div class=\"prompt output_prompt\">ZpFP" + no + "</div>\n"
            s += "        <div class=\"output_text output_subarea output_execute_result\">\n";
            s += "          <div class=\"timings\">%.3f&nbsp;s</div>\n" % TmeFP[n]
            s += "          <pre>" + htmlspecialchars(str(OutFP[n])) + "</pre>\n"
            s += "        </div>\n";
            s += "      </div>\n";
            s += "    </div>\n";
            s += "  </div>\n";
        if OutLC[n] is not None:
            s += "  <div class=\"output_wrapper\">\n"
            s += "    <div class=\"output\">\n"
            s += "      <div class=\"output_area\">\n"
            s += "        <div class=\"prompt output_prompt\">ZpLC" + no + "</div>\n"
            s += "        <div class=\"output_text output_subarea output_execute_result\">\n";
            s += "          <div class=\"timings\">%.3f&nbsp;s</div>\n" % TmeLC[n]
            s += "          <pre>" + htmlspecialchars(str(OutLC[n])) + "</pre>\n"
            s += "        </div>\n";
            s += "      </div>\n";
            s += "    </div>\n";
            s += "  </div>\n";
        if OutLF[n] is not None:
            s += "  <div class=\"output_wrapper\">\n"
            s += "    <div class=\"output\">\n"
            s += "      <div class=\"output_area\">\n"
            s += "        <div class=\"prompt output_prompt\">ZpLF" + no + "</div>\n"
            s += "        <div class=\"output_text output_subarea output_execute_result\">\n";
            s += "          <div class=\"timings\">%.3f&nbsp;s</div>\n" % TmeLF[n]
            s += "          <pre>" + htmlspecialchars(str(OutLF[n])) + "</pre>\n"
            s += "        </div>\n";
            s += "      </div>\n";
            s += "    </div>\n";
            s += "  </div>\n";

        s += "</div>\n";

    ffoot = io.open("foot.html")
    if ffoot:
        s += ffoot.read()
        ffoot.close()

    if fh:
        fh.write(s)
        fh.close()
    else:
        return s


def to_latex(filename, outfile=None):
    fh = None
    if outfile is not None:
        fh = io.open(outfile, "w")
    code = readfile(filename)
    OutCR, TmeCR = my_exec(code, ZpCR, QpCR)
    OutFP, TmeFP = my_exec(code, ZpFP, QpFP)
    OutLC, TmeLC = my_exec(code, ZpLC, QpLC)
    OutLF, TmeLF = my_exec(code, ZpLF, QpLF)

    nn = 0
    s = ""
    for n in range(len(code)):
        # Input
        if code[n][0] == "%% markdown\n":
            no = ""
            #s += "Markdown\n\n"
        else:
            nn += 1
            no = "~[" + str(nn) + "]:"
            s += "\\begin{tabular}{rl}\n"
            s += "\\cIn\n"
            for bloc in code[n]:
                for line in bloc.split("\n"):
                    if line != "": s += " & \\verb?" + line + "? \\\\\n"
            # s += pygments.highlight("".join(code[n]), lexer, latex_formatter) + "\n\n"

            # Output
            if OutCR[n] is not None:
                s += "\\cZpCR\n"
                ans = str(OutCR[n])
                for line in ans.split("\n"):
                    if line != "": s += " & \\verb?" + line + "? \\\\\n"
            if OutFP[n] is not None:
                s += "\\cZpFP\n"
                ans = str(OutFP[n])
                for line in ans.split("\n"):
                    if line != "": s += " & \\verb?" + line + "? \\\\\n"
            if OutLC[n] is not None:
                s += "\\cZpLC\n"
                ans = str(OutLC[n])
                for line in ans.split("\n"):
                    if line != "": s += " & \\verb?" + line + "? \\\\\n"
            if OutLF[n] is not None:
                s += "\\cZpLF\n"
                ans = str(OutLF[n])
                for line in ans.split("\n"):
                    if line != "": s += " & \\verb?" + line + "? \\\\\n"

            s += "\\end{tabular}\n\n"
            s += "\\bigskip\n\n"

    if fh:
        fh.write(unicode(s))
        fh.close()
    else:
        return s
