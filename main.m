clear;clc;

% type 1 (exponential: constant Q10)
Equation1 = fittype('a .* exp(b .* x)', ...
    'independent', 'x', 'dependent', 'y', ...
    'coefficients', {'a', 'b'});
initialParams1 = [0.59, 0.11];
opts1 = fitoptions('Method', 'NonlinearLeastSquares', 'StartPoint', initialParams1);
opts1.Lower = [0, 0];
% [fitresult, goodness] = fit(x, y,Equation1, 'StartPoint', initialParams1);

% type 2 (% Sigmoid / Logistic)
Equation2 = fittype('L ./ (1 + exp(-k .* (x - x0)))',...
    'independent', 'x', 'dependent', 'y', ...
    'coefficients', {'L', 'k', 'x0'});% upper % slope % centre
initialParams2 = [12, 0.2, 21];
opts2 = fitoptions('Method', 'NonlinearLeastSquares', 'StartPoint', initialParams2);
opts2.Lower = [0, 0, -Inf];

% type 3 (log-polynomial)
Equation3 = fittype('exp(a + b.* x + c .* x.^2)', ...
    'independent', 'x', 'dependent', 'y', ...
    'coefficients', {'a', 'b', 'c'});
initialParams3 = [0.32, 0.1012, -0.0005 ];
opts3 = fitoptions('Method', 'NonlinearLeastSquares', 'StartPoint', initialParams3);
opts3.Upper = [Inf, Inf,0];
% fitresult = fit(x, y, Equation3, opts);
  

% dataN = readtable('Data-Fluxnet_Eurosites.xlsx','Sheet','normal','Range','A1:I14601','ReadVariableNames',true);
dataH = readtable('Data-Fluxnet_Eurosites.xlsx','Sheet','heat','Range','A1:I14601','ReadVariableNames',true);
sz = [40 7];
varTypes = {'double','string','double','double','double','double','double'};
varNames = {'No','Site','I_hitype','I_hir2','II_hitype','II_hir2', 'hire_norm_area'};
temps = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);

for i = 1:max(dataH.No)
    sitei = dataH(dataH.No==i,:);
    temps{i,1} = i;
    temps{i,2} = string(cell2mat(sitei.SITE(1)));    
    idx_valid = isfinite(sitei.Re)&isfinite(sitei.Ta);
    re = sitei.Re(idx_valid);
    ta = sitei.Ta(idx_valid);
    idx0 = re>0;
    re = re(idx0);
    ta = ta(idx0);
    idx_Tmax = find(ta==max(ta));
    idx_I = 1:idx_Tmax; % warming branch is from the beginning of the year to Tmax
    idx_II = idx_Tmax:length(ta); % cooling branch is from Tmax to the end of the year

%   calculate HI
    pI = [ta(idx_I), re(idx_I)];
    taI = pI(:,1); reI = pI(:,2);
    pI_minT = floor(min(taI)); pI_maxT = ceil(max(taI));
    movingT0_I = []; movingER0_I = [];
    binsI = (pI_maxT-pI_minT)*2;
    for p = 1:binsI
        TI = pI_minT + (p-1);
        datapoints = find(taI <= (TI+3) & taI >= TI);
        movingT0_I(p) = mean(taI(datapoints));
        movingER0_I(p) = mean(reI(datapoints));
    end
    movingT0_I = movingT0_I';movingER0_I = movingER0_I';
    valid_indicesI = ~isnan(movingT0_I) & ~isnan(movingER0_I);
    movingT_I = movingT0_I(valid_indicesI);
    movingER_I = movingER0_I(valid_indicesI);
    
    pII = [ta(idx_II), re(idx_II)];
    taII = pII(:,1); reII = pII(:,2);
    pII_minT = floor(min(taII)); pII_maxT = ceil(max(taII));
    movingT0_II = []; movingER0_II = [];
    binsII = (pII_maxT-pII_minT)*2;
    for p = 1:binsII
        TII = pII_minT + (p-1);
        datapoints = find(taII <= (TII+3) & taII >= TII);
        movingT0_II(p) = mean(taII(datapoints));
        movingER0_II(p) = mean(reII(datapoints));
    end
    movingT0_II = movingT0_II';movingER0_II = movingER0_II';
    valid_indicesII = ~isnan(movingT0_II) & ~isnan(movingER0_II);
    movingT_II = movingT0_II(valid_indicesII);
    movingER_II = movingER0_II(valid_indicesII);
    
    
%     find the best fit
    [fitResult1I, goodness1I] = fit(movingT_I, movingER_I,Equation1, opts1);
    rSquare1I = goodness1I.rsquare;
    [fitResult2I, goodness2I] = fit(movingT_I, movingER_I,Equation2, opts2);
    rSquare2I = goodness2I.rsquare;
    [fitResult3I, goodness3I] = fit(movingT_I, movingER_I,Equation3, opts3);
    rSquare3I = goodness3I.rsquare;
    
    rSquare_I = [rSquare1I, rSquare2I, rSquare3I];
    [max_rSquare_I, max_rSquare_index_I] = max(rSquare_I);
    type_I = max_rSquare_index_I;
    if type_I ==1
        coeff_valuesI = coeffvalues(fitResult1I);
        equationI = @(x) coeff_valuesI(1) .* exp(coeff_valuesI(2) .* x);
        fitResultI = fitResult1I;
    elseif type_I == 2
        coeff_valuesI = coeffvalues(fitResult2I);
        equationI = @(x) coeff_valuesI(1) ./ (1 + exp(-coeff_valuesI(2) .* (x - coeff_valuesI(3))));
        fitResultI = fitResult2I;
    elseif type_I == 3
        coeff_valuesI = coeffvalues(fitResult3I);
        equationI = @(x) exp(coeff_valuesI(1) + coeff_valuesI(2).* x + coeff_valuesI(3) .* x.^2);
        fitResultI = fitResult3I;
    end
    temps{i,3} = type_I;
    temps{i,4} = max_rSquare_I;
    
    [fitResult1II, goodness1II] = fit(movingT_II, movingER_II,Equation1, opts1);
    rSquare1II = goodness1II.rsquare;
    [fitResult2II, goodness2II] = fit(movingT_II, movingER_II,Equation2, opts2);
    rSquare2II = goodness2II.rsquare;
    [fitResult3II, goodness3II] = fit(movingT_II, movingER_II,Equation3, opts3);
    rSquare3II = goodness3II.rsquare;
    
    rSquare_II = [rSquare1II, rSquare2II, rSquare3II];
    [max_rSquare_II, max_rSquare_index_II] = max(rSquare_II);
    type_II = max_rSquare_index_II;
    if type_II ==1
        coeff_valuesII = coeffvalues(fitResult1II);
        equationII = @(x) coeff_valuesII(1) .* exp(coeff_valuesII(2) .* x);
        fitResultII = fitResult1II;
    elseif type_II == 2
        coeff_valuesII = coeffvalues(fitResult2II);
        equationII = @(x) coeff_valuesII(1) ./ (1 + exp(-coeff_valuesII(2) .* (x - coeff_valuesII(3))));
        fitResultII = fitResult2II;
    elseif type_II == 3
        coeff_valuesII = coeffvalues(fitResult3II);
        equationII = @(x) exp(coeff_valuesII(1) + coeff_valuesII(2).* x + coeff_valuesII(3) .* x.^2);
        fitResultII = fitResult3II;
    end
    temps{i,5} = type_II;
    temps{i,6} = max_rSquare_II;
    
    delta_T = 0.1;
    x_values = floor(min(ta)):delta_T:ceil(max(ta));
    y_valuesI = feval(fitResultI, x_values);
    y_valuesII = feval(fitResultII, x_values);
    hys_area = sum(y_valuesI-y_valuesII)*delta_T;
    hys_area_norm = hys_area/(x_values(end)-x_values(1))/max(abs([y_valuesI; y_valuesII]));
    temps{i,7} = hys_area_norm;
    
end
writetable(temps,'HEAT_sites_HIRE.csv','WriteRowNames',true,'Delimiter',',','QuoteStrings',true)
