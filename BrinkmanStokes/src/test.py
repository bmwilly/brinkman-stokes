def bar(a, b):
    return a + b


def foo(a, b):
    c = bar(a, b)
    print(c)
    return c


d = 1
e = foo(d, 2)


foo(1, 2)
