import sys

def pw_char(i):
    if i % 4 == 0:
        return "-"
    elif i % 4 == 1:
        return "/"
    elif i % 4 == 2:
        return "|"
    else:
        return "\\"

for i in range(1000000):
    sys.stdout.write("{}\r".format(pw_char(i)))
