from sage.rings.polynomial.polynomial_ring import PolynomialRing_general
from sage.matrix.matrix_space import MatrixSpace
import io

import pygments
lexer = pygments.lexers.PythonLexer()
formatter = pygments.formatters.HtmlFormatter()

randints = [ ]
def random_element(ring, prec=None, degree=None):
    global random_seed
    random_seed += 1
    for _ in range(len(randints), random_seed+1):
        randints.append(ZZ.random_element(2^1000))
    if isinstance(ring, PolynomialRing_general):
        if degree is None:
            degree = 2
        R = ring.base_ring()
        coeffs = [ random_element(R, prec=prec) for _ in range(degree+1) ]
        return ring(coeffs).change_ring(R.fraction_field())
    elif isinstance(ring, MatrixSpace):
        size = ring.ncols() * ring.nrows()
        R = ring.base_ring()
        coeffs = [ random_element(R, prec=prec, degree=degree) for _ in range(size) ]
        return ring(coeffs)
    else:
        if prec is None:
            return ring(randints[random_seed])
        else:
            return ring(randints[random_seed]).add_bigoh(prec)


def my_exec(code, ring=Zp):
    global random_seed
    random_seed = -1;
    Zp = ring
    Out = [ ]
    for cmd in code:
        Out.append(None)
        length = len(cmd)
        if length == 0: continue
        cmd_exec = ""
        for n in range(length-1):
            cmd_exec += preparse(cmd[n])
        cmd_eval = preparse(cmd[-1])
        exec(cmd_exec)
        try:
            Out[-1] = eval(cmd_eval)
        except SyntaxError:
            exec(cmd_eval)
    return Out
    

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
    OutCR = my_exec(code, ZpCR)
    OutFP = my_exec(code, ZpFP)
    OutLP = my_exec(code, ZpLP)
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
        if OutLP[n] is not None:
            s += "ZpLP:\n" + str(OutLP[n]) + "\n"
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
    OutCR = my_exec(code, ZpCR)
    OutFP = my_exec(code, ZpFP)
    OutLP = my_exec(code, ZpLP)

    s = ""
    for n in range(len(code)):
        no = "&nbsp;[" + str(n+1) + "]:"
        s += " <div class=\"cell border-box-sizing code_cell rendered\">\n"

        # Input
        s += "  <div class=\"input\">\n"
        s += "    <div class=\"prompt input_prompt\">In" + no + "</div>\n"
        s += "    <div class=\"inner_cell\">\n"
        s += "      <div class=\"input_area\">\n"
        s += pygments.highlight("".join(code[n]), lexer, formatter) + "\n"
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
            s += "          <pre>" + htmlspecialchars(str(OutFP[n])) + "</pre>\n"
            s += "        </div>\n";
            s += "      </div>\n";
            s += "    </div>\n";
            s += "  </div>\n";
        if OutLP[n] is not None:
            s += "  <div class=\"output_wrapper\">\n"
            s += "    <div class=\"output\">\n"
            s += "      <div class=\"output_area\">\n"
            s += "        <div class=\"prompt output_prompt\">ZpLP" + no + "</div>\n"
            s += "        <div class=\"output_text output_subarea output_execute_result\">\n";
            s += "          <pre>" + htmlspecialchars(str(OutLP[n])) + "</pre>\n"
            s += "        </div>\n";
            s += "      </div>\n";
            s += "    </div>\n";
            s += "  </div>\n";

        s += "</div>\n";

    if fh:
        fh.write(s)
        fh.close()
    return s
