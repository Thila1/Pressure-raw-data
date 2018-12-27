clc;
clear all;
close all;

%============ Consider note for getFitdataGraph function ==================================
%           'poly1','poly2' 'poly3','poly4','poly5','poly6','poly7','poly8','poly9' <==need
% ploy_type :  1       2       3       4       5       6       7       8       9
%
 ploy_type=1; % <==ploynamial type , is defind above <====================== check this pls
%==========================================================================================


%======================Note of Read Ms Excel Files==============================
% [timeBub1,vol_bub1,press_bub1,timeBub2,vol_bub2,press_bub2,frequncy,photo,uv]=allData_Pressure_time_volume_voltage_frequncy(url)
%  URL refers to the link of the file (file location to be defined)
% './' to consider the current folder as the root folder
% The files to be analyzed is referred through the folder
%===============================================================================
%--------The file names to be analyzed
%  0_5MNacl_0ml.xlsx   
%  0_5MNacl_2ml.xlsx
%  0_5MNacl_4ml.xlsx
%  0_5MNacl_6ml.xlsx
%  2MNacl_0ml.xlsx
%  2MNacl_2ml.xlsx  
%  2MNacl_4ml.xlsx
%  2MNacl_6ml.xlsx
%  distill_0ml.xlsx
%  distill_2ml.xlsx
%  distill_4ml.xlsx
%  distill_6ml.xlsx


%Make the necessary changes accordingly  (plot graphs use below ,between of below function definds and this line)
% Change the URL of the file to analyse other data 

nameOfChemical='Distilled water';%<=== fill this data (automaticaly add to title of graph)---------
M_value='';% concentration of the solution
ml_amount='6ml';
analysis_is_distill=1;  %If distilled water please change it to 1. if not leave it at 0.
file_url='./Analysis/distill_6ml.xlsx'; % ./ current working folder is root

is_orginal_only_transis_frequncy_curve=1;    %<= orginal data point graph only no need curve fit equation then 
is_orginal_only_uv_frequncy_curve=1;           % is_orginal_only_transis_frequncy_curve=1;
                                             % Or need curve fit graph
                                             % is_orginal_only_transis_frequncy_curve=0;

[timeBub1,vol_bub1,press_bub1,timeBub2,vol_bub2,press_bub2,frequncy,photo_Voltage,uv_Voltage]=allData_Pressure_time_volume_voltage_frequncy(file_url,analysis_is_distill);
%----------------------------------------------------------------------------------------
percentage='(97% ENA)';
percentage2='(with 97% ENA)';
if(analysis_is_distill==1)
    percentage='';
    percentage2='';
end


figure, % figure 1 Pressure (Pa) vs Time (s)
plot(timeBub1,press_bub1,'o','color','red');
hold on
plot(timeBub2,press_bub2,'o','color','blue');
grid on;
title({['Pressure (Pa) vs Time (s) for the cavitation bubble with ',ml_amount,' air drawn out for'], [M_value,' ',nameOfChemical,' solution ',percentage2]});
ylabel('Pressure (Pa)');
xlabel('Time (s)');
legend('Bubble one','Bubble two','Orientation','horizontal');
hold off;


figure, % figure 2 Pressure (Pa) vs Volume (mm3)
plot(vol_bub1,press_bub1,'c*','color','red');
hold on;
plot(vol_bub2,press_bub2,'c*','color','blue');
grid on;
title({['Pressure (Pa) vs Volume (mm3) for the cavitation bubble with ',ml_amount,' air drawn out for'], [M_value,' ',nameOfChemical,' solution ',percentage2]});
ylabel('Pressure (Pa)');
xlabel('Volume (mm3)');
legend(['Bubble one ',ml_amount],['Bubble two ',ml_amount],'Orientation','horizontal');
hold off;


[fit_trans,goodness_data_tran]=getFitdataGraph(frequncy,photo_Voltage); % fitting function
% figure 3 process
figure, % figure 3 Voltage (mV) vs Frequency (kHz) photo Transistor
if (goodness_data_tran.rsquare>0.8700)&&(is_orginal_only_transis_frequncy_curve==0)
    transis_data=plot(fit_trans,frequncy,photo_Voltage);
else
    plot(frequncy,photo_Voltage);
end
grid on;
title({['Voltage (mV) vs Frequency (kHz) for ',M_value],[' ',nameOfChemical,' ',percentage,' with ',ml_amount,' air drawn out for the phototransistor']});
ylabel('Frequency (kHz)');
xlabel('Voltage (mV)');
legend('Photo Transistor','Location','northwest');
hold off;
% 

[fit_uv,goodness_data_uv]=getFitdataGraph(frequncy,uv_Voltage); % fitting function
% figure 4 process
figure, % figure 4 Voltage (mV) vs Frequency (kHz) UV Diod
if (goodness_data_uv.rsquare>0.8700)&&(is_orginal_only_uv_frequncy_curve==0)
    uv_data=plot(fit_uv,frequncy,uv_Voltage);
else
    plot(frequncy,uv_Voltage);
end
grid on;
title({['Voltage (mV) vs Frequency (kHz) for ',M_value],[' ',nameOfChemical,' ',percentage,' with ',ml_amount,' air drawn out for the UV doide']});
ylabel('Frequency (kHz)');
xlabel('Voltage (mV)');
legend('UV diode','Location','northwest');
hold off;
% 
% 
listofLegents=[];
figure, % figure 5 Voltage (mV) vs Frequency for both
if(goodness_data_tran.rsquare>0.8700)&&(is_orginal_only_transis_frequncy_curve==0)
    plot(transis_data(2).XData,transis_data(2).YData);
else
    plot(frequncy,photo_Voltage);%fit_trans
end
  hold on;
if(goodness_data_uv.rsquare>0.8700)&&(is_orginal_only_uv_frequncy_curve==0)
    plot(uv_data(2).XData,uv_data(2).YData);  
else
    plot(frequncy,uv_Voltage);
end
grid on;
title({['Voltage (mV) vs Frequency (kHz) for ', M_value],[' ',nameOfChemical,' ',percentage,' with ',ml_amount,' air drawn out']});
ylabel('Frequency (kHz)');
xlabel('Voltage (mV)');
legend('Photo Transistor','UV diode','Location','northwest');
hold off;

% display data
if(is_orginal_only_transis_frequncy_curve==0)
%     if (goodness_data_tran.sse<1)
        disp('Transis graph Coffiee');
        disp(fit_trans);
        disp(goodness_data_tran);
%     else
%         disp('Not match Equation to Best fitting');
%         disp(goodness_data_tran);
%     end
end
if(is_orginal_only_uv_frequncy_curve==0)
%     if (goodness_data_uv.sse<1)
        disp('UV graph Coffiee');
        disp(fit_uv);
        disp(goodness_data_uv);
%     else
%         disp('Not match Equation to Best fitting');
%         disp(goodness_data_uv);
%     end
else
    disp('Not Curve of Equation fit Orginal data point graph only');
end



%=========== Below function definds =========================================

%-----create fit curve function using ploynomial type---------------------------
function [fit_data,goodness_data]=getFitdataGraph(x_axis_data,y_axis_mean) % this function is create fit value of seleted ploynomial type
%     poly_set={'poly1','poly2','poly3','poly4','poly5','poly6','poly7','poly8','poly9'};% ploy nomial type
    fo = fitoptions('Method','NonlinearLeastSquares','Lower',[0,0],'Upper',[Inf,max(x_axis_data)]);
    ft = fittype('(a*(x^2) + b*x + c) /((x^3) + e*(x^2) + f*x+d)','options',fo);
    [fit_data,goodness_data]= fit(x_axis_data',y_axis_mean',ft,'normalize','on');
    
%     ploy=poly_set{ploy_type}; 'Normalize','on'
end
%----------------create fit function end ---------------------------------------

%----------------------data read functions in excels-----------------------
function [timeBub1,vol_bub1,press_bub1,timeBub2,vol_bub2,press_bub2,frequncy,photo_Voltage,uv_Voltage]=allData_Pressure_time_volume_voltage_frequncy(url,analysis_distill) % this function implemt given below read_data.. two functions
    timeBub1=[]; % exsis bubble time for bub 1
    timeBub2=[]; % exsis bubble time for bub 2
    
    vol_bub1=[]; % exsis bubble volum for bub 1
    vol_bub2=[]; % exsis bubble volum for bub 2
    
    press_bub1=[]; % exsis bubble pressure for bub 1
    press_bub2=[]; % exsis bubble pressure for bub 2
    
    sheet1_name='sheet1';
    sheet2_name='sheet2';
    if(analysis_distill==1)
         sheet1_name='Bubble analysis';
         sheet2_name='sono analysis';
    end
    
    if(analysis_distill==0)
          sheet1_name='sheet1';
          sheet2_name='sheet2';
    end
    
    % read time , volum , pressure used below defind function
    [witch_row_start,timeR,vol_bub1R,vol_bub2R,press_bub1R,press_bub2R]=read_data_Pressure_time_volume(url,sheet1_name);% read sheet1 datas
    
    % read frequncy , reading of photo transistor ,reading of UV diod pressure used below defind function
    [frequncy,photo_Voltage,uv_Voltage]=read_data_voltage_frequncy(url,sheet2_name);% read sheet 2 datas
    
    % removing NaN & null values
    for r=witch_row_start:1:length(timeR)
        if ~(vol_bub1R(r)==0)
            if(~isnan(vol_bub1R(r)))&&(~isnan(press_bub1R(r)))&&(~isnan(timeR(r)))
                timeBub1=[timeBub1,timeR(r)];
                vol_bub1=[vol_bub1,vol_bub1R(r)];
                press_bub1=[press_bub1,press_bub1R(r)];
            end
        end
        if ~(vol_bub2R(r)==0)
             if(~isnan(vol_bub2R(r)))&&(~isnan(press_bub2R(r)))&&(~isnan(timeR(r)))
                timeBub2=[timeBub2,timeR(r)];
                vol_bub2=[vol_bub2,vol_bub2R(r)];
                press_bub2=[press_bub2,press_bub2R(r)];
             end
        end
    end
    %
end
function [witch_row_start,timeR,vol_bub1R,vol_bub2R,press_bub1R,press_bub2R]=read_data_Pressure_time_volume(url,sheet1_name)% read sheet 1 datas
    [num,text,~]= xlsread(url,sheet1_name); % read excel file

    [~,col]=size(text); % read number of colums and rows

    vol_bub1_col=0; % colum number of vol bub1
    vol_bub2_col=0; % colum number of vol bub2
    press_bub1_col=0; % colum number of press bub1
    press_bub2_col=0; % colum number of press bub2
    time_col=0;% colum number of time

    state=0;% set all colum number state
    witch_row_start=0;
    for j=4:1:10
        for i=1:1:col
            if(isequal(char(cellstr(cell(text(j,i)))),'Vol bub 1 V1 (mm3)'))
                vol_bub1_col=i;
                state=state+1;
                witch_row_start=j;
            end
            if(isequal(char(cellstr(cell(text(j,i)))),'Vol bub 2 V2 (mm3)'))
                vol_bub2_col=i;
                state=state+1;
            end
            if(isequal(char(cellstr(cell(text(j,i)))),'Pressure in bubble 1'))
                press_bub1_col=i;
                state=state+1;
            end
            if(isequal(char(cellstr(cell(text(j,i)))),'Pressure in bubble 2'))
                press_bub2_col=i;
                state=state+1;
            end
            if(isequal(char(cellstr(cell(text(j,i)))),'Time'))
                time_col=i;
                state=state+1;
            end
        end
        if(state==5)
            break;
        end
    end
    if(state==5)
        press_bub1R=[];
        press_bub2R=[];
        
        timeR=num(:,time_col);
        vol_bub1R=num(:,vol_bub1_col);
        vol_bub2R=num(:,vol_bub2_col);
        
        if(press_bub1_col>0)&&(~((max(vol_bub1R))==0))
            press_bub1R=num(:,press_bub1_col);
        end
        if(press_bub2_col>0)&&(~((max(vol_bub2R))==0))
            press_bub2R=num(:,press_bub2_col);
        end
    end
end
function [frequncy,photo_Voltage,uv_Voltage]=read_data_voltage_frequncy(url,sheet2_name)% read sheet 2 datas
    [num,text,~]= xlsread(url,sheet2_name);
    
    [rows,col]=size(text); % read number of colums and rows

    photoTran_col=0; % colum number of photo transistor input
    uvDiod_col=0; % colum number of uv diod input
    frequncy_col=0; % colum number of frequncy

    states=0;% set all colum number state
    witch_row=0;
    for j=4:1:rows
        for i=1:1:col
            if(isequal(char(cellstr(cell(text(j,i)))),'Phototransistor Reading (mV)'))
                photoTran_col=i;
                states=states+0;
                witch_row=j;
            end
            if(isequal(char(cellstr(cell(text(j,i)))),'UV diode  (mV)'))
                uvDiod_col=i;
                states=states+1;
            end
            if(isequal(char(cellstr(cell(text(j,i)))),'Frequency (kHz)'))
                frequncy_col=i;
                states=states+1;
            end
            
        end
        if(states==3)
            break;
        end
    end
    uv_Voltage=[];
    photo_Voltage=[];
    frequncy=[];
     if(states==3)
        % remove nan values
        for r=witch_row:1:length(num(:,frequncy_col))
            if(~isnan(num(r,uvDiod_col)))
                uv_Voltage=[uv_Voltage,num(r,uvDiod_col)];
            end
            if(~isnan(num(r,photoTran_col)))
                photo_Voltage=[photo_Voltage,num(r,photoTran_col)];
            end
            if(~isnan(num(r,frequncy_col)))
                frequncy=[frequncy,num(r,frequncy_col)];
            end
        end
     end
end
%--------------------------------------------------------------------------