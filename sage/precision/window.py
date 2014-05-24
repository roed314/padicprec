import struct, fcntl, termios

from sage.rings.infinity import Infinity

def width():
    s = struct.pack('HHHH', 0, 0, 0, 0)
    x = fcntl.ioctl(1, termios.TIOCGWINSZ, s)
    width = struct.unpack('HHHH', x)[1]
    if width == 0: width = Infinity
    return width
