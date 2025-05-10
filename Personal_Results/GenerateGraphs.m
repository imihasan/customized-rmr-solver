%% Generate graphs used in the paper
% This script generates graphs used in the paper
%
% Author: Ibrahim Mohammed I. Hasan (imihasan@kth), 2025

%% Write text files
% navigate the path to where this file is saved
path=pwd;

folders=dir(fullfile(path,"S8R","With Scaling"));
for f=1:length(folders)-2
    trials=dir(fullfile(path,"S8R", "With Scaling",string(folders(f+2).name)+"\",'*.mat'));
    for t=1:length(trials)
        load(fullfile(path,"S8R", "With Scaling",string(folders(f+2).name),trials(t).name));
        l=0:length(xsol)-1;
        time=(1/frequency_solution)*l';
        vars=[time xsol(:,1:33) xsol(:,34:end) force_vector norm_fv_in_ground norm_fv_rotated norm_Vec_H2GCC Vec_GC2GEE Vec_H2GCC];
        varsNames=["time", muscle_order(1:33), muscle_order(34:end), "Fx", "Fy", "Fz", "Fgx", "Fgy", "Fgz", "Fnx", "Fny", "Fnz", "HCx", "HCy", "HCz", "GC2GEx", "GC2GEy", "GC2GEz", "H2GCCx", "H2GCCy", "H2GCCz"];
        Tab = array2table(vars,'VariableNames',varsNames);
        writetable(Tab, fullfile(path, "S8R", "With Scaling",string(folders(f+2).name),trials(t).name(1:end-4) +".txt"));
        fprintf("Writing file: "+ string(folders(f+2).name) + " - " +trials(t).name(1:end-4) +"\n")
    end
end
%% Plot contact force magnitude
close all
path=fullfile(pwd,"S8R", "With Scaling");

models="senstivity_planar"; %Change here which stability models you want to plot according to the conditions in the blelow if statement

if models == "constraints"
    folders=["Point", "Circle", "Polynomial","Ellipse"];
    Conditions=["Orthoload" folders];
    color=[0.6350 0.0780 0.1840; 0 0.4470 0.7410; 0.4660 0.6740 0.1880;  0.9290 0.6940 0.1250 ];
    trials_names=["Lateral Raise", "Posterior Raise", "Anterior Raise"];

elseif models== "penalties"
    folders=["Conditional", "Planar", "Curve","None"];
    Conditions=["Orthoload" folders];
    color=[0.8500 0.3250 0.0980; 0.3010 0.7450 0.9330; 0.4940 0.1840 0.5560; 1 0 1];
    trials_names=["Lateral Raise", "Posterior Raise", "Anterior Raise"];

elseif models=="senstivity_curve"
    folders=["Curve w=1", "Curve", "Curve w=10"];
    Conditions=["Orthoload" , "Looser", "Normal", "Tighter"];
    trials_names=["Lateral Raise", "Posterior Raise", "Anterior Raise"];
    color=[0.4660 0.6740 0.1880; 0 0.4470 0.7410; 0.6350 0.0780 0.1840];

elseif models=="senstivity_conditional"
    folders=["Conditional w=1", "Conditional", "Conditional w=10"];
    Conditions=["Orthoload" , "Looser", "Normal", "Tighter"];
    trials_names=["Lateral Raise"];
    color=[0.4660 0.6740 0.1880; 0 0.4470 0.7410; 0.6350 0.0780 0.1840];

elseif models=="senstivity_planar"
    folders=["Planar w=1", "Planar", "Planar w=10"];
    Conditions=["Orthoload" , "Looser", "Normal", "Tighter"];
    trials_names=["Lateral Raise","Posterior Raise","Anterior Raise"];
    color=[0.4660 0.6740 0.1880; 0 0.4470 0.7410; 0.6350 0.0780 0.1840];

elseif models=="senstivity_polynomial"
    folders=["Polynomial Loose", "Polynomial", "Polynomial Tight"];
    Conditions=["Orthoload" , "Looser", "Normal", "Tighter"];
    trials_names=["Lateral Raise","Posterior Raise","Anterior Raise"];
    color=[0.4660 0.6740 0.1880; 0 0.4470 0.7410; 0.6350 0.0780 0.1840];

elseif models=="all"
    folders=["Point", "Circle", "Polynomial","Ellipse", "Conditional", "Planar", "Curve","None"];
    Conditions=["Orthoload" folders];
    color=[0.6350 0.0780 0.1840; 0 0.4470 0.7410; 0.4660 0.6740 0.1880;  0.9290 0.6940 0.1250; 0.8500 0.3250 0.0980; 0.3010 0.7450 0.9330; 0.4940 0.1840 0.5560; 1 0 1 ];
    trials_names=["Lateral Raise", "Posterior Raise", "Anterior Raise"];

elseif models=="rotator1"
    folders=["Point", "Planar"];
    Conditions=["Orthoload" folders];
    color=[0.6350 0.0780 0.1840; 0 0.4470 0.7410; 0.4940 0.1840 0.5560; 1 0 1 ];
    trials_names=["Lateral Raise", "Posterior Raise", "Anterior Raise"];

elseif models=="rotator2"
    folders=["Planar", "Curve"];
    Conditions=["Orthoload" folders];
    color=[ 0.3010 0.7450 0.9330; 0.4940 0.1840 0.5560];
    trials_names=["Lateral Raise", "Posterior Raise", "Anterior Raise"];

elseif models=="scaled"
    folders=["Planar", "Test"];
    Conditions=["Orthoload" folders];
    color=[ 0.3010 0.7450 0.9330; 0.4940 0.1840 0.5560];
    trials_names=["Lateral Raise"];
end

trials=dir(fullfile(path,folders(1)+"\",'*.mat'));



fs=2000;
fc=200;
[b,a]=butter(6,fc/(fs/2));

for t= 1:length(trials)
    %Read Orthoload
    Orthoload=readtable(fullfile(path, Conditions(1),trials(t).name(1:end-4) +".txt"));
    OrtholoadG=readtable(fullfile(path, Conditions(1),"F Glenoid",trials(t).name(1:end-4) +".txt"));
    
    %Find Sync indices
    indx1=find(Orthoload.Marker==1,1);
    indx2=find(Orthoload.Marker(indx1:end)==0,1)-1;

    t_ortho=Orthoload.Time(indx1:indx2)-Orthoload.Time(indx1);
    mag_ortho=filtfilt(b,a,sqrt(Orthoload.Fx(indx1:indx2).^2+Orthoload.Fy(indx1:indx2).^2+Orthoload.Fz(indx1:indx2).^2));
    fig=figure;
    plot(t_ortho,mag_ortho*100/830,'k',LineWidth=3)
    hold on
    for c=1:length(folders)

        cond=readtable(fullfile(path, folders(c),trials(t).name(1:end-4) +".txt"));

        if trials_names(t)=="Posterior Raise"
            idx2_osim=find(cond.time>=3.4,1);
        elseif trials_names(t)=="Lateral Raise"
            idx2_osim=find(cond.time>=4,1);
        else
            idx2_osim=find(cond.time>=4.5,1);
        end
        force=filtfilt(b,a,sqrt(cond.Fx.^2+cond.Fy.^2+cond.Fz.^2));
        force=force(1:idx2_osim);
        plot(cond.time(1:idx2_osim),force*100/830,LineWidth=3,Color=[color(c,:)])
        error{t,c}=(force-spline(t_ortho,mag_ortho,cond.time(1:idx2_osim)))*100/830;
        RMSE{t,c}=rmse(force,spline(t_ortho,mag_ortho,cond.time(1:idx2_osim)))*100/830;

        rotator{t,c}=[(cond.Subscapularis_S(1:idx2_osim) + cond.Subscapularis_M(1:idx2_osim) + cond.Subscapularis_I(1:idx2_osim))/3, ...
            (cond.Infraspinatus_I(1:idx2_osim)+cond.Infraspinatus_S(1:idx2_osim))/2, (cond.Supraspinatus_P(1:idx2_osim)+cond.Supraspinatus_A(1:idx2_osim))/2,...
            cond.TeresMinor(1:idx2_osim)];
    end
    ylabel("GH-JCF [% BW]",FontSize=15)
    xlabel("Time (s)",FontSize=15)
    ax = gca;
    ax.FontSize = 15; 
    box off
    title(trials_names(t),FontSize=20)
    % legend(Conditions,Location="best",FontSize=20)
    if trials_names(t)=="Posterior Raise"
        xlim([0 3.4])
        ylim([0 400])
    elseif trials_names(t)=="Lateral Raise"
        xlim([0 4])
        ylim([0 400])
    else 
        xlim([0 4.5])
        
    end
    ylim([0 400])
    saveas(fig,fullfile(pwd,"Figures Scaled",trials_names(t)+"_"+models+"_GHJCFMagnitude.pdf")); %name your file here
end

%% RPlot root mean square error (RMSE)
trials_namess=["L Raise", "P Raise", "A Raise"];
for t=1:length(trials_names)
    for c=1:length(Conditions)-1
        Max(t,c)=max(error{t,c});
        Min(t,c)=min(error{t,c});
        RRMSE(t,c)=RMSE{t,c};
    end 
end

fig=figure;
Cat= categorical(trials_namess);
Cat=reordercats(Cat,trials_namess);
hb=bar(Cat,Max,'BarWidth',0.3);

hold on
hb1=bar(Cat,Min,'BarWidth',0.3);

for k = 1:numel(hb)  
    xtips = hb(k).XEndPoints;
    ytips = hb(k).YEndPoints;
    hb(k).FaceColor=color(k,:);
    hb(k).EdgeColor=color(k,:);
    hb(k).LineWidth=2;
    hold on
end

for k = 1:numel(hb1)  
    xtips = hb1(k).XEndPoints;
    ytips = hb1(k).YEndPoints;
    hb1(k).FaceColor=color(k,:);
    hb1(k).EdgeColor=color(k,:);
    hb1(k).LineWidth=2;
    hold on
end

for k = 1:numel(hb)                                                      
    xtips = hb(k).XEndPoints;
    ytips = hb(k).YEndPoints;
    scatter(xtips,RRMSE(:,k)',30,color(k,:),'o','filled','MarkerEdgeColor','k',LineWidth=2);
end

%legend(Conditions(2:end),Location="best",FontSize=20)
ylabel("Estimation Error [% BW]")
ylim([-200 300])
ax = gca;
ax.FontSize = 15;
box off
saveas(fig,fullfile(pwd,"Figures Scaled",trials_names(t)+"_"+models+"_GHJCFMagnitude.RMSE.pdf"));

%% Rotator Cuff Activations - Timeseries

for t=1:length(trials_names)
    cond=readtable(fullfile(path, folders(c),trials(t).name(1:end-4) +".txt"));
    if trials_names(t)=="Posterior Raise"
        idx2_osim=find(cond.time>=3.4,1);
    elseif trials_names(t)=="Lateral Raise"
        idx2_osim=find(cond.time>=4,1);
    else
        idx2_osim=find(cond.time>=4.5,1);
    end
    fig=figure;
    hold on
    for c=1:length(Conditions)-1
        plot(cond.time(1:idx2_osim),rotator{t,c}(:,1),color=color(c,:),LineWidth=2,LineStyle="-")
        plot(cond.time(1:idx2_osim),rotator{t,c}(:,2),color=color(c,:),LineWidth=2,LineStyle=":")
        plot(cond.time(1:idx2_osim),rotator{t,c}(:,3),color=color(c,:),LineWidth=2,LineStyle="--")
        % plot(cond.time(1:idx2_osim),rotator{t,c}(:,4),color=color(c,:),LineWidth=2,LineStyle="-.")
    end 
    %legend('Subscapularis','Infraspinatus','Supraspinatus','Teres Minor')
    box off
    ax = gca;
    ax.FontSize = 20;
    title(trials_names(t))
    ylabel("Activation")
    xlabel("Time (s)")
    ylim([0 1])
    if trials_names(t)=="Posterior Raise"
        xlim([0 3.4])
    elseif trials_names(t)=="Lateral Raise"
        xlim([0 4])
    else 
        xlim([0 4.5])
        
    end
    saveas(fig,fullfile(pwd,"Figures Scaled",trials_names(t)+"_"+models+"_Rotator_Cuff_Activation.pdf"));
end

%% Rotator Cuff Activation - Bar Plot

fig=figure;
trials_namess=["L Raise", "P Raise", "A Raise"];
Cat= categorical(trials_namess);
Cat=reordercats(Cat,trials_namess);
hb=bar(Cat,m_rotator,'EdgeColor','k','BarWidth',0.4);


color1 = [0.9 0.7 1 ; 0.8 1 1;
   0.9 1 0.1; 0.9 0.8 0.3];

hold on

for k = 1:numel(hb)  
    xtips = hb(k).XEndPoints;
    ytips = hb(k).YEndPoints;
    hb(k).FaceColor=color(k,:);
    hb(k).EdgeColor='k';

    hold on
    for i=1:length(xtips)
        for s=1:4
            bl=scatter(xtips(i),m_rotator_indiv{i,k}(s),100,color1(s,:),'filled','o','MarkerEdgeColor','k',LineWidth=1);
        end
        
    end
end
 
ylabel("Mean Activation",fontsize=15)
box off
ax = gca;
ax.FontSize = 15;
h = get(gca,'Children');

ylim([0 0.5])
title(models)
saveas(fig,fullfile(pwd,"Figures Scaled",models+"_RotatorCuffActivation.pdf")); %name your file here

%% GH-JCF Direction
matlab.graphics.internal.setPrintPreferences('DefaultPaperPositionMode','auto')
set(groot,'defaultFigurePaperPositionMode','auto')

models="senstivity_polynomial"; %Change here 

if models == "constraints"
    folders=["Point", "Circle", "Polynomial","Ellipse"];
    Conditions=["Orthoload" folders];
    color=[0.6350 0.0780 0.1840; 0 0.4470 0.7410; 0.4660 0.6740 0.1880;  0.9290 0.6940 0.1250 ];
    trials_names=["Lateral Raise", "Posterior Raise", "Anterior Raise"];
elseif models== "penalties"
    folders=["Conditional", "Planar", "Curve"];
    Conditions=["Orthoload" folders];
    color=[0.8500 0.3250 0.0980; 0.3010 0.7450 0.9330; 0.4940 0.1840 0.5560; 1 0 1];
    trials_names=["Lateral Raise", "Posterior Raise", "Anterior Raise"];
elseif models=="senstivity_curve"
    folders=["Curve w=1", "Curve", "Curve w=10"];
    Conditions=["Orthoload" , "Looser", "Normal", "Tighter"];
    trials_names=["Lateral Raise", "Posterior Raise", "Anterior Raise"];
    color=[0.4660 0.6740 0.1880; 0 0.4470 0.7410; 0.6350 0.0780 0.1840];
elseif models=="senstivity_conditional"
    folders=["Conditional w=1", "Conditional", "Conditional w=10"];
    Conditions=["Orthoload" , "Looser", "Normal", "Tighter"];
    trials_names=["Lateral Raise"];
    color=[0.4660 0.6740 0.1880; 0 0.4470 0.7410; 0.6350 0.0780 0.1840];
elseif models=="senstivity_planar"
    folders=["Planar w=1", "Planar", "Planar w=10"];
    Conditions=["Orthoload" , "Looser", "Normal", "Tighter"];
    trials_names=["Anterior Raise", "Anterior Raise", "Anterior Raise"];
    color=[0.4660 0.6740 0.1880; 0 0.4470 0.7410; 0.6350 0.0780 0.1840];

elseif models=="senstivity_polynomial"
    folders=["Polynomial Loose", "Polynomial", "Polynomial Tight"];
    Conditions=["Orthoload" , "Looser", "Normal", "Tighter"];
    trials_names=["Lateral Raise","Posterior Raise","Anterior Raise"];
    color=[0.4660 0.6740 0.1880; 0 0.4470 0.7410; 0.6350 0.0780 0.1840];
elseif models=="all"
    folders=["Point", "Circle", "Polynomial","Ellipse", "Conditional", "Planar", "Curve","None"];
    Conditions=["Orthoload" folders];
    color=[0.6350 0.0780 0.1840; 0 0.4470 0.7410; 0.4660 0.6740 0.1880;  0.9290 0.6940 0.1250; 0.8500 0.3250 0.0980; 0.3010 0.7450 0.9330; 0.4940 0.1840 0.5560; 1 0 1 ];
    trials_names=["Lateral Raise", "Posterior Raise", "Anterior Raise"];
end



trials=dir(fullfile(path,folders(1)+"\",'*.mat'));



for t=1:length(trials_names)
    OrtholoadG=readtable(fullfile(path, Conditions(1),"F Glenoid",trials(t).name(1:end-4) +".txt"));
    Orthoload=readtable(fullfile(path, Conditions(1),trials(t).name(1:end-4) +".txt"));
    
    indx1=find(Orthoload.Marker==1,1);
    if trials_names(t)=="Posterior Raise"
        indx2=find(Orthoload.Time(indx1:end)-Orthoload.Time(indx1)>=3.4,1)-1+indx1;
    elseif trials_names(t)=="Lateral Raise"
        indx2=find(Orthoload.Time(indx1:end)-Orthoload.Time(indx1)>=4,1)-1+indx1;
    else 
        indx2=find(Orthoload.Time(indx1:end)-Orthoload.Time(indx1)>=4.5,1)-1+indx1;
    end

    load(fullfile(path, folders(1),trials(t).name(1:end-4) +".mat"));
    t_ortho=Orthoload.Time(indx1:indx2)-Orthoload.Time(indx1);

    F_ortho=[OrtholoadG.Fgx(indx1:indx2) OrtholoadG.Fgy(indx1:indx2) OrtholoadG.Fgz(indx1:indx2)];
    for f=1:length(F_ortho)
        cosThetaa = max(min(dot([0 0 1],F_ortho(f,:))/(norm([0 0 1])*norm(F_ortho(f,:))),1),-1);

        F_ortho_norm(f,:)=(25/abs(cosThetaa))*F_ortho(f,:)/(norm(F_ortho(f,:)));
    end
    fig=figure;
   
    theta=0:0.001:2*pi;
    sc=norm(Vec_H2GCC(end,:))*1000;

    %Ellipse
    x=sc*0.34*cos(theta);
    y=sc*0.61*sin(theta);
    
    %Polynomial
    theta=rad2deg(theta);
    L=(-5.430378773797311e-12*theta.^6+5.880441464447251e-09*theta.^5+-2.309728860559176e-06*theta.^4+3.899916118022262e-04*theta.^3+-0.024526120335986*theta.^2+0.116047310713347*theta+50.915928515941545)/100;
    theta=deg2rad(theta);
    xl=sc*L.*cos(theta);
    yl=sc*L.*sin(theta);
    
    %uncomment if you want to plot the ellipse and the polynomial
    %perimeters
    %plot(x,y,LineWidth=1,LineStyle="--",Color=[color(4,:) 1]) % Plot Ellipse Perimeter
    
    hold on
    plot(2*yl,2*xl,LineWidth=1,LineStyle="--",Color=color(1,:)) % Plot Polynomial Fit
    plot(yl,xl,LineWidth=1,LineStyle="--",Color=color(2,:)) % Plot Polynomial Fit
    plot(0.5*yl,0.5*xl,LineWidth=1,LineStyle="--",Color=color(3,:)) % Plot Polynomial Fit

    axis equal
    ax=gca;
    hold on
    scatter(F_ortho_norm(:,1),F_ortho_norm(:,2), 40, [0.5 0.5 0.5], 'filled','MarkerFaceAlpha',0.2)
    scatter(mean(F_ortho_norm(:,1)),mean(F_ortho_norm(:,2)), 150, [0.5 0.5 0.5], 'filled','LineWidth',1,'MarkerEdgeColor','k')
    
    for c=1:length(Conditions)-1
        load(fullfile(path, folders(c),trials(t).name(1:end-4) +".mat"));
        t_optim(t,c)=tOptim;

        l=0:length(xsol)-1;
        time=(1/frequency_solution)*l';

        %Crop variables to the correct time instant
        if trials_names(t)=="Posterior Raise"
            idx2_osim=find(time==3.4);
        elseif trials_names(t)=="Lateral Raise"
            idx2_osim=find(time==4);
        else
            idx2_osim=find(time==4.5);
        end


        for i=1:length(norm_Vec_H2GCC)
            xhat=-norm_Vec_H2GCC(i,:);
            norm_Vec_GC2GEE=Vec_GC2GEE(i,:)/norm(Vec_GC2GEE(i,:));
            Vxdot=dot(norm_Vec_GC2GEE,xhat);
            Vy=norm_Vec_GC2GEE-Vxdot*xhat;
            yhat=Vy/norm(Vy);
            zhat=cross(xhat,yhat);
        
            R=[xhat; yhat; zhat];
            cosTheta = dot(Vec_H2GCC(i,:),force_vector(i,:))/(norm(Vec_H2GCC(i,:))*norm(force_vector(i,:)));
            norm_fv_rotated(i,:) = R*norm_fv_in_ground(i,:)'*(norm(1000*Vec_H2GCC(i,:))/cosTheta);

        end
        
        scatter(-norm_fv_rotated(1:idx2_osim,2), -norm_fv_rotated(1:idx2_osim,3), 40, color(c,:), 'filled','MarkerFaceAlpha',0.2);
        scatter(mean(-norm_fv_rotated(1:idx2_osim,2)), mean(-norm_fv_rotated(1:idx2_osim,3)), 150, color(c,:), 'filled','LineWidth',1,'MarkerEdgeColor','k');
        

        % difference_rmse_y(t,c)=rmse(-norm_fv_rotated(1:idx2_osim,3),spline(t_ortho,F_ortho_norm(:,2),time));
        % 
        % difference_rmse_x(t,c)=rmse(-norm_fv_rotated(1:idx2_osim,2),spline(t_ortho,F_ortho_norm(:,1),time));

        clear norm_fv_rotated force_vector
        
    end
    clear F_ortho_norm F_ortho
        vv=R*Vec_GC2GEE(end,:)';
        h2gc=R*Vec_H2GCC(end,:)'*1000;
        radius = norm(h2gc)*0.5;

        %norm(Vec_GC2GEE(end,:))/norm(Vec_H2GCC(end,:))
        p=nsidedpoly(1000, 'Center', [0,0], 'Radius', radius);
        
        %plot(p, 'FaceColor', 'none',LineWidth=1,EdgeColor=[0 0.4470 0.7410],EdgeAlpha=1,LineStyle='--')
        % axis equal
        box off
        % leg(1:2:2*length(Conditions(2:end)))=Conditions(2:end);
        % legend(["" "Experimental" "" "None" "" "Estimated"],'Location','best')

        h2ge=R*Vec_GC2GEE(end,:)'*1000;        
        scatter(-h2gc(2), -h2gc(3),100,'x',MarkerEdgeColor='k',LineWidth=2)

        h = gca;
        xlabel("X [mm]")   % corresponding roughly to OpenSim X axis (horizontal pointing forward)
        ylabel("Y [mm]")   % corresponding to OpenSim Y axis (vertical pointing upwards)
        ax = gca;
        ax.FontSize = 15; 
        xlim([-35 35])
        ylim([-50 50])
        title(trials_names(t))
        saveas(fig,fullfile(pwd,"Figures Scaled",trials_names(t)+"_"+models+"_GHJCF_Direction.pdf")); %name the figure here
end

%% RMSE of direction
trials_namess=["L Raise", "P Raise", "A Raise"];

fig=figure;
Cat= categorical(trials_namess);
Cat=reordercats(Cat,trials_namess);
hb=bar(Cat,(difference_rmse_y+difference_rmse_x)/2,'BarWidth',0.7);


for k = 1:numel(hb)  
    xtips = hb(k).XEndPoints;
    ytips = hb(k).YEndPoints;
    hb(k).FaceColor=color(k,:);
    hb(k).EdgeColor=color(k,:);
    %hb(k).LineWidth=2;

end
 
ylabel("RMSE [mm]",fontsize=15)
box off
ax = gca;
ax.FontSize = 20;
h = get(gca,'Children');

saveas(fig,fullfile(pwd,"Figures Scaled","_"+models+"_RMSE_Direction.pdf")); %name the figure here

%% Computational time

trials_namess=["Lateral Raise", "P Raise"];

fig=figure;
Cat= categorical(trials_namess);
Cat=reordercats(Cat,trials_namess);
hb=bar(Cat,t_optim(1:2,1:8),'BarWidth',0.7);


for k = 1:numel(hb)  
    xtips = hb(k).XEndPoints;
    ytips = hb(k).YEndPoints;
    hb(k).FaceColor=color(k,:);
    hb(k).EdgeColor=color(k,:);
    %hb(k).LineWidth=2;

end

ylabel("Computational Time [s]",fontsize=15)
box off
ax = gca;
ax.FontSize = 20;
h = get(gca,'Children');

saveas(fig,fullfile(pwd,"Figures Scaled","_"+models+"_Computational_Time.pdf")); %name the figure here

%% GH-JCF Direction with no Stability Condition (None)
folders=["None"];
trials=dir(fullfile(path,folders(1)+"\",'*.mat'));
trials_names=["Lateral Raise","Posterior Raise","Anterior Raise"];
Conditions=["Orthoload", "None"];
color=[13, 74, 33; 120, 0, 26; 98, 152, 210; 0, 0, 97]/255;

fig=figure;
for t=1:length(trials_names)

    OrtholoadG=readtable(fullfile(path, Conditions(1),"F Glenoid",trials(t).name(1:end-4) +".txt"));
    Orthoload=readtable(fullfile(path, Conditions(1),trials(t).name(1:end-4) +".txt"));
    
    indx1=find(Orthoload.Marker==1,1);
    if trials_names(t)=="Posterior Raise"
        indx2=find(Orthoload.Time(indx1:end)-Orthoload.Time(indx1)>=3.4,1)-1+indx1;
    elseif trials_names(t)=="Lateral Raise"
        indx2=find(Orthoload.Time(indx1:end)-Orthoload.Time(indx1)>=4,1)-1+indx1;
    else 
        indx2=find(Orthoload.Time(indx1:end)-Orthoload.Time(indx1)>=4.5,1)-1+indx1;
    end
    

    F_ortho=[OrtholoadG.Fgx(indx1:indx2) OrtholoadG.Fgy(indx1:indx2) OrtholoadG.Fgz(indx1:indx2)];
    for f=1:length(F_ortho)
        F_ortho_norm(f,:)=25*F_ortho(f,:)/norm(F_ortho(f,:));
    end

    % figure
    hold on
    for c=1:length(Conditions)-1
        load(fullfile(path, folders(c),trials(t).name(1:end-4) +".mat"));


        for i=1:length(norm_Vec_H2GCC)
            xhat=-norm_Vec_H2GCC(i,:);
            norm_Vec_GC2GEE=Vec_GC2GEE(i,:)/norm(Vec_GC2GEE(i,:));
            Vxdot=dot(norm_Vec_GC2GEE,xhat);
            Vy=norm_Vec_GC2GEE-Vxdot*xhat;
            yhat=Vy/norm(Vy);
            zhat=cross(xhat,yhat);
        
            R=[xhat; yhat; zhat];
        
            cosTheta = dot(Vec_H2GCC(i,:),force_vector(i,:))/(norm(Vec_H2GCC(i,:))*norm(force_vector(i,:)));
            norm_fv_rotated(i,:) = R*norm_fv_in_ground(i,:)'*(norm(1000*Vec_H2GCC(i,:))/cosTheta);
        
        end
        scatter(-norm_fv_rotated(:,2), -norm_fv_rotated(:,3), 40, color(t,:), 'filled','MarkerFaceAlpha',0.1);
        hold on
        scatter(mean(-norm_fv_rotated(:,2)), mean(-norm_fv_rotated(:,3)), 150, color(t,:), 'filled','LineWidth',1,'MarkerEdgeColor','k');
        
        
    end 

        axis equal
        box off
        leg(1:2:2*length(Conditions(2:end)))=Conditions(2:end);
        
        h = gca;
        xlabel("X [mm]")   % corresponding roughly to OpenSim X axis (horizontal pointing forward)
        ylabel("Y [mm]")   % corresponding to OpenSim Y axis (vertical pointing upwards)
        ax = gca;
        ax.FontSize = 20; 

end
vv=R*Vec_GC2GEE(end,:)';
radius = norm([vv(2) vv(3)])*1000;
p=nsidedpoly(1000, 'Center', [0,0], 'Radius', radius);
p1=nsidedpoly(1000, 'Center', [0,0], 'Radius', 25);
plot(p, 'FaceColor', 'none',LineWidth=1,LineStyle='--',EdgeColor=[0 0.4470 0.7410])
title("None Model")
ax=gca;
ax.FontSize=20;
saveas(fig,fullfile(pwd,"Figures Scaled","Loc_none.pdf")); %name the figure here
axis equal

%% GH-JCF Direction Sensitivity

close all
models="constraints";

if models=="penalties"
    folders=["Curve w=1", "Curve", "Curve w=10"];
    Conditions=["Orthoload" , "Looser", "Normal", "Tighter"];
elseif models == "constraints"
    folders=["Polynomial Tight", "Polynomial", "Polynomial Loose"];
    Conditions=["Orthoload" , "Looser", "Normal", "Tighter"];
end


trials=dir(fullfile(path,folders(1)+"\",'*.mat'));
trials_names=["Lateral Raise"];


color=[0.4660 0.6740 0.1880; 0 0.4470 0.7410; 0.6350 0.0780 0.1840;   0.9290 0.6940 0.1250 ;  0.8500 0.3250 0.0980; 0.3010 0.7450 0.9330; 0.4940 0.1840 0.5560; 1 0 1];

for t=1:length(trials_names)
    OrtholoadG=readtable(fullfile(path, Conditions(1),"F Glenoid",trials(t).name(1:end-4) +".txt"));
    Orthoload=readtable(fullfile(path, Conditions(1),trials(t).name(1:end-4) +".txt"));
    
    indx1=find(Orthoload.Marker==1,1);
    if trials_names(t)=="Posterior Raise"
        indx2=find(Orthoload.Time(indx1:end)-Orthoload.Time(indx1)>=3.4,1)-1;
    elseif trials_names(t)=="Lateral Raise"
        indx2=find(Orthoload.Time(indx1:end)-Orthoload.Time(indx1)>=4,1)-1;
    else 
        indx2=find(Orthoload.Time(indx1:end)-Orthoload.Time(indx1)>=4.5,1)-1;
    end

    load(fullfile(path, folders(1),trials(t).name(1:end-4) +".mat"));

    F_ortho=[OrtholoadG.Fgx(indx1:indx2) OrtholoadG.Fgy(indx1:indx2) OrtholoadG.Fgz(indx1:indx2)];
    for f=1:length(F_ortho)
        cosTheta = max(min(dot([0 0 -1],F_ortho(f,:))/(norm([0 0 -1])*norm(F_ortho(f,:))),1),-1);
        F_ortho_norm(f,:)=25*F_ortho(f,:)/(norm(F_ortho(f,:)));
    end
    fig=figure;
    theta=0:0.001:2*pi;

    L=(-5.430378773797311e-12*theta.^6+5.880441464447251e-09*theta.^5+-2.309728860559176e-06*theta.^4+3.899916118022262e-04*theta.^3+-0.024526120335986*theta.^2+0.116047310713347*theta+50.915928515941545)/100;
    
    
    sc=norm(Vec_H2GCC(end,:))*1000;

    xc=sc*sqrt(0.5)*cos(theta);
    yc=sc*sqrt(0.5)*sin(theta);
    
    x=sc*0.34*cos(theta);
    y=sc*0.61*sin(theta);
    
    theta=rad2deg(theta);
    L=(-5.430378773797311e-12*theta.^6+5.880441464447251e-09*theta.^5+-2.309728860559176e-06*theta.^4+3.899916118022262e-04*theta.^3+-0.024526120335986*theta.^2+0.116047310713347*theta+50.915928515941545)/100;
    
    theta=deg2rad(theta);
    xl=sc*L.*cos(theta);
    yl=sc*L.*sin(theta);
    
    hold on
    plot(yl,xl,LineWidth=2,LineStyle=":",Color=[color(2,:) 1]) %uncomment
    % if you want to draw the polynomial perimeter
    plot(2*yl,2*xl,LineWidth=2,LineStyle=":",Color=[color(1,:) 1])
    plot(0.5*yl,0.5*xl,LineWidth=2,LineStyle=":",Color=[color(3,:) 1])

    axis equal
    ax=gca;
    hold on
    scatter(F_ortho_norm(:,1),F_ortho_norm(:,2), 40, [0.5 0.5 0.5], 'filled','MarkerFaceAlpha',0.2)
    scatter(mean(F_ortho_norm(:,1)),mean(F_ortho_norm(:,2)), 150, [0.5 0.5 0.5], 'filled','LineWidth',1,'MarkerEdgeColor','k')
    clear F_ortho_norm F_ortho
    for c=1:length(Conditions)-1
        load(fullfile(path, folders(c),trials(t).name(1:end-4) +".mat"));
        l=0:length(xsol)-1;
        time=(1/frequency_solution)*l';

        %Crop variables to the correct time instant
        if trials_names(t)=="Posterior Raise"
            idx2_osim=find(time==3.4);
        elseif trials_names(t)=="Lateral Raise"
            idx2_osim=find(time==4);
        else
            idx2_osim=find(time==4.5);
        end


        for i=1:length(norm_Vec_H2GCC)
            xhat=-norm_Vec_H2GCC(i,:);
            norm_Vec_GC2GEE=Vec_GC2GEE(i,:)/norm(Vec_GC2GEE(i,:));
            Vxdot=dot(norm_Vec_GC2GEE,xhat);
            Vy=norm_Vec_GC2GEE-Vxdot*xhat;
            yhat=Vy/norm(Vy);
            zhat=cross(xhat,yhat);
        
            R=[xhat; yhat; zhat];
        
            cosTheta = dot(Vec_H2GCC(i,:),force_vector(i,:))/(norm(Vec_H2GCC(i,:))*norm(force_vector(i,:)));
            norm_fv_rotated(i,:) = R*norm_fv_in_ground(i,:)'*(norm(1000*Vec_H2GCC(i,:))/cosTheta);
       
        end

        scatter(-norm_fv_rotated(1:idx2_osim,2), -norm_fv_rotated(1:idx2_osim,3), 40, color(c,:), 'filled','MarkerFaceAlpha',0.2);
        scatter(mean(-norm_fv_rotated(1:idx2_osim,2)), mean(-norm_fv_rotated(1:idx2_osim,3)), 150, color(c,:), 'filled','LineWidth',1,'MarkerEdgeColor','k');
        clear norm_fv_rotated force_vector 
    end 
        vv=R*Vec_GC2GEE(end,:)';
        radius = norm([vv(2) vv(3)])*1000;


        p=nsidedpoly(1000, 'Center', [0,0], 'Radius', radius);
        %plot(p, 'FaceColor', 'none',LineWidth=1,EdgeColor=[0 0.4470 0.7410],EdgeAlpha=1,LineStyle='--')
        axis equal
        box off
        leg(1:2:2*length(Conditions(2:end)))=Conditions(2:end);

        h2gc=R*Vec_H2GCC(end,:)'*1000;
        scatter(-h2gc(2), -h2gc(3),100,'x',MarkerEdgeColor='k',LineWidth=2)

        h = gca;
        xlabel("X [mm]")   % corresponding roughly to OpenSim X axis (horizontal pointing forward)
        ylabel("Y [mm]")      % corresponding to OpenSim Y axis (vertical pointing upwards)
        ax = gca;
        ax.FontSize = 18; 
        ylim([-50 50])
        xlim([-35 35])
        title(trials_names(t))
        saveas(fig,fullfile(pwd,"Figures Scaled",trials_names(t)+"_Sensitivity_"+models+"NT.pdf")); %name the figure here
end

%% Plot Reserve Actuators vs Joint Moments
matlab.graphics.internal.setPrintPreferences('DefaultPaperPositionMode','auto')
set(groot,'defaultFigurePaperPositionMode','auto')
close all
path="C:\Users\ibrah\OneDrive - KTH\Model Validation - Orthoload\rmr-solver-1.0 - original\rmr-solver-1.0\Personal_Results\S8R\With GH Constraint\New Transformation";

models="all"; %Change here 

if models == "constraints"
    folders=["Point", "Circle", "Polynomial","Ellipse"];
    Conditions=["Orthoload" folders];
    color=[0.6350 0.0780 0.1840; 0 0.4470 0.7410; 0.4660 0.6740 0.1880;  0.9290 0.6940 0.1250 ];
    trials_names=["Lateral Raise", "Posterior Raise", "Anterior Raise"];
elseif models== "penalties"
    folders=["Conditional", "Planar", "Curve","None"];
    Conditions=["Orthoload" folders];
    color=[0.8500 0.3250 0.0980; 0.3010 0.7450 0.9330; 0.4940 0.1840 0.5560; 1 0 1];
    trials_names=["Lateral Raise", "Posterior Raise", "Anterior Raise"];

elseif models=="all"
    folders=["None","Point", "Circle", "Polynomial","Ellipse", "Conditional", "Planar", "Curve"];
    Conditions=["Orthoload" folders];
    color=[1 0 1 ; 0.6350 0.0780 0.1840; 0 0.4470 0.7410; 0.4660 0.6740 0.1880;  0.9290 0.6940 0.1250; 0.8500 0.3250 0.0980; 0.3010 0.7450 0.9330; 0.4940 0.1840 0.5560];
    trials_names=["Lateral Raise", "Posterior Raise", "Anterior Raise"];

end

trials=dir(fullfile(path,folders(1)+"\",'*.mat'));
trial=["MSabdkg1", "MSbfkg1", "MSffkg1"];

coordinate=["scapula_abduction", "scapula_elevation", "scapula_upward_rot", ...
    "scapula_winging", "plane_elv", "shoulder_elv", "axial_rot", "elbow_flexion"];

coordinatess=["Scapula Abduction", "Shoulder Elevation", "Scapula Upward Rot", ...
    "Scapula Winging", "Plane of Elevation", "Shoulder Elevation", "Humeral Rot", "Elbow flex"];

fs=2000;
fc=200;
[b,a]=butter(6,fc/(fs/2));

for t= 1:length(trials)
    r_act=importdata("C:\Users\ibrah\OneDrive - KTH\Model Validation - Orthoload\rmr-solver-1.0 - original\rmr-solver-1.0\Personal_Results\S8R\With GH Constraint\New Transformation\Curve\inverse_dynamics_"+trial(t)+".sto");
    clo_head=string(r_act.colheaders)';
    
    fig=figure;
    % fig.WindowState="maximized";
    for cord=1:length(coordinate)
        act=r_act.data(:,find(clo_head==coordinate(cord)+"_moment"));
        
        subplot(2,4,cord)
        plot(r_act.data(:,1),act,'k',LineWidth=2,Color=[0.5 0.5 0.5])
        hold on

        for c=1:length(folders)
            cond=readtable(fullfile(path, folders(c),trials(t).name(1:end-4) +".txt"));
    
            if trials_names(t)=="Posterior Raise"
                idx2_osim=find(cond.time>=3.4,1);
                xlim([0 3.4])
            elseif trials_names(t)=="Lateral Raise"
                idx2_osim=find(cond.time>=4,1);
                 xlim([0 4])
            else
                idx2_osim=find(cond.time>=4.5,1);
                xlim([0 4.5])
            end
  
            act_osim=eval(['cond.' char(coordinate(cord))]);
            plot(cond.time(1:idx2_osim),act_osim(1:idx2_osim),LineWidth=2,Color=[color(c,:)])
            ylabel("Moment [N.m]",FontSize=8)
            xlabel("Time (s)",FontSize=8)
            ax = gca;
            ax.FontSize =8; 
            box off
            title(coordinatess(cord),FontSize=8,Interpreter="none")
            %legend("Joint", "R Actuator")
        end
    end
mkdir(fullfile(pwd,"\Figures Scaled\ReserveActuators"));
saveas(fig,fullfile(pwd,"Figures Scaled\ReserveActuators",trials_names(t)+"_"+models+"_ReserveActuator.pdf"));

end