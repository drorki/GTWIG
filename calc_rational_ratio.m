function rational_ratios = calc_rational_ratio(hMax)

val = zeros(hMax^2, 1);
nn = 0;
for ii = 1:hMax
    for jj = 1:hMax
        nn = nn+1;
        val(nn) = ii/jj;
        %disp([num2str(ii) '/' num2str(jj) ' = ' num2str(val(nn))])
    end
end

rational_ratios = unique(val);
rational_ratios( rational_ratios==1 ) = [];