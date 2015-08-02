msizes = [4,5,6,7,8]; mt = zeros(length(msizes),1);
i = 1;
for msize = msizes
    etoc = mv_fun();
    fprintf('vectors generated in %8.3e seconds\n',etoc)
    mt(i) = etoc;
    i = i + 1;
end

