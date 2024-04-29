close all;clear;clc;

for pic_id=1:1:26
    load('data.mat','I1');
    list=I1(927,1177:1192);
    xx=1:1:length(list);
    
    [fitresult, gof] = createFit(xx, list);
    calibration_result(pic_id)=fitresult.c;
end
% calibration_result=unwrap2(calibration_result);
xf=0:10:250;
yf=calibration_result-calibration_result(1);
figure,scatter(xf,yf)

delete_index=[5,7,15,17,25];xf(delete_index)=[];yf(delete_index)=[];
yf=unwrap2(yf);
figure,scatter(xf,yf)


function [fitresult, gof] = createFit(xx, list)
%CREATEFIT(XX,LIST)

[xData, yData] = prepareCurveData( xx, list );

ft = fittype( 'a*sin(b*x+c)+d', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
a0=(max(list)-min(list))/2;d0=mean(list);
b0=2*pi/10;c0=asin((list(1)-d0)/a0);
opts.StartPoint = [a0,b0,c0,d0];
opts.Lower=[a0*0.1,b0*0.1,-pi,d0*0.1];
opts.Upper=[a0*10 ,b0*10 , pi,d0*10 ];

[fitresult, gof] = fit( xData, yData, ft, opts );

figure( 'Name', '拟合' );
h = plot( fitresult, xData, yData );
legend( h, 'list vs. xx', '拟合', 'Location', 'NorthEast', 'Interpreter', 'none' );
xlabel( 'xx', 'Interpreter', 'none' );
ylabel( 'list', 'Interpreter', 'none' );
% ylim([2.2,3.4])
grid on
end
