function bar(a, b)
    return a + b
end

function foo(a, b)
    c = bar(a, b)
    println(c)
end

d = 1
e = foo(d, 2)


foo(1, 2)
