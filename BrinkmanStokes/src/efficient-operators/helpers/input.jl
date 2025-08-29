function user_input(prompt::String = "")
    print(prompt)
    return chomp(readline())
end
