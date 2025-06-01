%Sasha Kramer
%MBARI

%%%Linear relationship between POC flux and specific phytoplankton taxa
%Figures 4-7 in Kramer et al., 2025

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%INITAL DATA PROCESSING
cd ~/18S/QIIME/Jan23
%load DNA data: taxonomy, metadata, feature table
%Isolate samples to use

%Load taxonomy and separate out phototrophs
cd ~/18S/taxonomy
load taxonomy.mat

for i = 1:17835
    k = taxo.Kingdom(i)==photo.Kingdom;
    s = taxo.Supergroup(i)==photo.Supergroup;
    d = taxo.Division(i)==photo.Division;
    c = taxo.Class(i)==photo.Class;
    ksdc = k+s+d+c;
    check = ksdc==4;
    if sum(check)>0
        indP(i,1) = 1;
    else indP(i,1) = 0;
    end
    clear k s d c ksdc check
end
clear i

%Only phototrophs
taxo_p = taxo(indP==1,:);
feat_p = feat(indP==1,:);

%Separate groups of samples by type and remove ASVs with no reads
NP_meta = meta([1:348,799:917,973:1026],:);
NP_featP = feat_p(:,[1:348,799:917,973:1026]);
NP_taxoP = taxo_p;

NP_taxoP(sum(NP_featP,2)==0,:) = [];
NP_featP(sum(NP_featP,2)==0,:) = [];

NA_meta = meta([349:798,918:972,1027:1080],:);
NA_featP = feat_p(:,[349:798,918:972,1027:1080]);
NA_taxoP = taxo_p;

NA_taxoP(sum(NA_featP,2)==0,:) = [];
NA_featP(sum(NA_featP,2)==0,:) = [];

clear feat feat_p indP meta photo taxo taxo_p

%Remove ASVs detected <20 reads max in any one sample
NP_taxoP(max(NP_featP,[],2)<20,:) = [];
NP_featP(max(NP_featP,[],2)<20,:) = [];

NA_taxoP(max(NA_featP,[],2)<20,:) = [];
NA_featP(max(NA_featP,[],2)<20,:) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Figure 4 analysis and generate figure
%Isolate ASVs from specific groups of samples to compare composition
%Surface NP
surf_NP = NP_featP(:,466:519);
surftax_NP = NP_taxoP;
surftax_NP(sum(surf_NP,2)==0,:) = [];
surf_NP(sum(surf_NP,2)==0,:) = [];
surfsum_NP = sum(surf_NP,2);

for i = 1:526
    if surftax_NP.Genus(i)=='NaN' && surftax_NP.Class(i)=='NaN'
        surft_NP(i,1) = 'Other';
    elseif surftax_NP.Genus(i)=='NaN'
        surft_NP(i,1) = surftax_NP.Class(i);
    else surft_NP(i,1) = surftax_NP.Genus(i);
    end
end
clear i

surftax_NP.Class(surftax_NP.Class=='NaN') = 'Other';
[~,indT] = sort(surftax_NP.Class);
surftax_NP = surftax_NP(indT,:);
surft_NP = surft_NP(indT,:);
surfsum_NP = surfsum_NP(indT,:);
surf_NP = surf_NP(indT,:);
clear indT

%Repeat for NP bulk and NP particles

%Find ASVs that are shared between the surface + bulk or surface + particles
[a1,~] = ismember(surftax_NP.FeatureID,bulktax_NP.FeatureID);
[a2,~] = ismember(surftax_NP.FeatureID,parttax_NP.FeatureID);
NPsurf_smpart = categorical(526,1);
for i = 1:526
    if a1(i)==0 && a2(i)==0
        NPsurf_smpart(i,1) = 'None';
    elseif a1(i)==1 && a2(i)==1
        NPsurf_smpart(i,1) = 'large';
    elseif a1(i)==0 && a2(i)==1
        NPsurf_smpart(i,1) = 'large';
    elseif a1(i)==1 && a2(i)==0
        NPsurf_smpart(i,1) = 'small';
    end
end
clear a1 a2

[NPsurf_taxo,indA] = sortrows(surftax_NP,[4,5]);
NPsurf_asvs = surf_NP(indA,:);
NPsurf_smpartS = NPsurf_smpart(indA,:);
clear indA

%Divide into 6 pigment-based groups:
%NPsurf_frac = summed ASVs in each group

%Surface NA
surf_NA = NA_featP(:,506:559);
surftax_NA = NA_taxoP;
surftax_NA(sum(surf_NA,2)==0,:) = [];
surf_NA(sum(surf_NA,2)==0,:) = [];
surfsum_NA = sum(surf_NA,2);

for i = 1:404
    if surftax_NA.Genus(i)=='NaN' && surftax_NA.Class(i)=='NaN'
        surft_NA(i,1) = 'Other';
    elseif surftax_NA.Genus(i)=='NaN'
        surft_NA(i,1) = surftax_NA.Class(i);
    else surft_NA(i,1) = surftax_NA.Genus(i);
    end
end
clear i

surftax_NA.Class(surftax_NA.Class=='NaN') = 'Other';
[~,indT] = sort(surftax_NA.Class);
surftax_NA = surftax_NA(indT,:);
surft_NA = surft_NA(indT,:);
surfsum_NA = surfsum_NA(indT,:);
surf_NA = surf_NA(indT,:);
clear indT

%Repeat for NA bulk and NA particles

%Find ASVs that are shared between the surface + bulk or surface + particles
[b1,~] = ismember(surftax_NA.FeatureID,bulktax_NA.FeatureID);
[b2,~] = ismember(surftax_NA.FeatureID,parttax_NA.FeatureID);
NAsurf_smpart = categorical(404,1);
for i = 1:404
    if b1(i)==0 && b2(i)==0
        NAsurf_smpart(i,1) = 'None';
    elseif b1(i)==1 && b2(i)==1
        NAsurf_smpart(i,1) = 'large';
    elseif b1(i)==0 && b2(i)==1
        NAsurf_smpart(i,1) = 'large';
    elseif b1(i)==1 && b2(i)==0
        NAsurf_smpart(i,1) = 'small';
    end
end
clear b1 b2

[NAsurf_taxo,indA] = sortrows(surftax_NA,[4,5]);
NAsurf_asvs = surf_NA(indA,:);
NAsurf_smpartS = NAsurf_smpart(indA,:);
clear indA

%Divide into 6 pigment-based groups:
%NAsurf_frac = summed ASVs in each group

%Bar plots of diversity in surface vs. small/large particles
%Separate into 6 pigment-based phytoplankton groups
%Look at ASVs in large, small, or None particles
%NP_smpbar, NA_smpbar

NP_smpbar_asv_frac = NP_smpbar_asv./sum(NP_smpbar_asv,2);
NA_smpbar_asv_frac = NA_smpbar_asv./sum(NA_smpbar_asv,2);

NP_smpbar_frac = NP_smpbar./sum(NP_smpbar,2);
NA_smpbar_frac = NA_smpbar./sum(NA_smpbar,2);

%Define colors to use for light and dark based on treemaps
colb = hex2rgb('#337539');
colb(2,:) = hex2rgb('#b8e0bb');
colb(3,:) = [1 1 1];

cold = hex2rgb('#5da899');
cold(2,:) = hex2rgb('#cee5e0');
cold(3,:) = [1 1 1];

colh = hex2rgb('#94caec');
colh(2,:) = hex2rgb('#dfeff9');
colh(3,:) = [1 1 1];

colg = hex2rgb('#c26a77');
colg(2,:) = hex2rgb('#edd2d6');
colg(3,:) = [1 1 1];

cols = hex2rgb('#9f4a96');
cols(2,:) = hex2rgb('#e5c6e1');
cols(3,:) = [1 1 1];

colo = hex2rgb('#7e2954');
colo(2,:) = hex2rgb('#e6b1cc');
colo(3,:) = [1 1 1];

%Plot results for each basin (diversity and relative seq abun)
figure(4),clf
subplot(2,2,1)
h = barh(1,NP_smpbar_frac(1,:),'stacked','EdgeColor',colb(1,:),'LineWidth',1.5);
for i = 1:3
    set(h(i),'FaceColor',colb(i,:))
end
clear i h
hold on
h = barh(2,NP_smpbar_frac(2,:),'stacked','EdgeColor',cold(1,:),'LineWidth',1.5);
for i = 1:3
    set(h(i),'FaceColor',cold(i,:))
end
clear i h
h = barh(3,NP_smpbar_frac(3,:),'stacked','EdgeColor',colh(1,:),'LineWidth',1.5);
for i = 1:3
    set(h(i),'FaceColor',colh(i,:))
end
clear i h
h = barh(4,NP_smpbar_frac(4,:),'stacked','EdgeColor',colg(1,:),'LineWidth',1.5);
for i = 1:3
    set(h(i),'FaceColor',colg(i,:))
end
clear i h
h = barh(5,NP_smpbar_frac(5,:),'stacked','EdgeColor',cols(1,:),'LineWidth',1.5);
for i = 1:3
    set(h(i),'FaceColor',cols(i,:))
end
clear i h
h = barh(6,NP_smpbar_frac(6,:),'stacked','EdgeColor',colo(1,:),'LineWidth',1.5);
for i = 1:3
    set(h(i),'FaceColor',colo(i,:))
end
clear i h
text([1.02,1.02,1.02,1.02,1.02,1.02],[1,2,3,4,5,6],{'55','276','88','37','29','39'},'fontsize',20)
axis ij
yticks([1:1:6])
yticklabels({'Diatoms','Dinoflagellates','Hacrobia','Chlorophytes','Dictyos+Pelagos','Other Ochrophyta'})
set(gca,'fontsize',20)
xlabel('Proportion of exported ASV diversity in surface samples')
axis([0 1 0 7])
title('North Pacific')

subplot(2,2,2)
h = barh(1,NA_smpbar_frac(1,:),'stacked','EdgeColor',colb(1,:),'LineWidth',1.5);
for i = 1:3
    set(h(i),'FaceColor',colb(i,:))
end
clear i h
hold on
h = barh(2,NA_smpbar_frac(2,:),'stacked','EdgeColor',cold(1,:),'LineWidth',1.5);
for i = 1:3
    set(h(i),'FaceColor',cold(i,:))
end
clear i h
h = barh(3,NA_smpbar_frac(3,:),'stacked','EdgeColor',colh(1,:),'LineWidth',1.5);
for i = 1:3
    set(h(i),'FaceColor',colh(i,:))
end
clear i h
h = barh(4,NA_smpbar_frac(4,:),'stacked','EdgeColor',colg(1,:),'LineWidth',1.5);
for i = 1:3
    set(h(i),'FaceColor',colg(i,:))
end
clear i h
h = barh(5,NA_smpbar_frac(5,:),'stacked','EdgeColor',cols(1,:),'LineWidth',1.5);
for i = 1:3
    set(h(i),'FaceColor',cols(i,:))
end
clear i h
h = barh(6,NA_smpbar_frac(6,:),'stacked','EdgeColor',colo(1,:),'LineWidth',1.5);
for i = 1:3
    set(h(i),'FaceColor',colo(i,:))
end
clear i h
text([1.02,1.02,1.02,1.02,1.02,1.02],[1,2,3,4,5,6],{'94','191','62','14','26','17'},'fontsize',20)
axis ij
yticks([1:1:6])
yticklabels({'Diatoms','Dinoflagellates','Hacrobia','Chlorophytes','Dictyos+Pelagos','Other Ochrophyta'})
set(gca,'fontsize',20)
xlabel('Proportion of exported ASV diversity in surface samples')
axis([0 1 0 7])
title('North Atlantic')

subplot(2,2,3)
h = barh(1,NP_smpbar_asv_frac(1,:),'stacked','EdgeColor',colb(1,:),'LineWidth',1.5);
for i = 1:3
    set(h(i),'FaceColor',colb(i,:))
end
clear i h
hold on
h = barh(2,NP_smpbar_asv_frac(2,:),'stacked','EdgeColor',cold(1,:),'LineWidth',1.5);
for i = 1:3
    set(h(i),'FaceColor',cold(i,:))
end
clear i h
h = barh(3,NP_smpbar_asv_frac(3,:),'stacked','EdgeColor',colh(1,:),'LineWidth',1.5);
for i = 1:3
    set(h(i),'FaceColor',colh(i,:))
end
clear i h
h = barh(4,NP_smpbar_asv_frac(4,:),'stacked','EdgeColor',colg(1,:),'LineWidth',1.5);
for i = 1:3
    set(h(i),'FaceColor',colg(i,:))
end
clear i h
h = barh(5,NP_smpbar_asv_frac(5,:),'stacked','EdgeColor',cols(1,:),'LineWidth',1.5);
for i = 1:3
    set(h(i),'FaceColor',cols(i,:))
end
clear i h
h = barh(6,NP_smpbar_asv_frac(6,:),'stacked','EdgeColor',colo(1,:),'LineWidth',1.5);
for i = 1:3
    set(h(i),'FaceColor',colo(i,:))
end
clear i h
axis ij
yticks([1:1:6])
yticklabels({'Diatoms','Dinoflagellates','Hacrobia','Chlorophytes','Dictyos+Pelagos','Other Ochrophyta'})
set(gca,'fontsize',20)
xlabel('Proportion of exported relative sequence abundance in surface samples')
axis([0 1 0 7])
title('North Pacific')

subplot(2,2,4)
h = barh(1,NA_smpbar_asv_frac(1,:),'stacked','EdgeColor',colb(1,:),'LineWidth',1.5);
for i = 1:3
    set(h(i),'FaceColor',colb(i,:))
end
clear i h
hold on
h = barh(2,NA_smpbar_asv_frac(2,:),'stacked','EdgeColor',cold(1,:),'LineWidth',1.5);
for i = 1:3
    set(h(i),'FaceColor',cold(i,:))
end
clear i h
h = barh(3,NA_smpbar_asv_frac(3,:),'stacked','EdgeColor',colh(1,:),'LineWidth',1.5);
for i = 1:3
    set(h(i),'FaceColor',colh(i,:))
end
clear i h
h = barh(4,NA_smpbar_asv_frac(4,:),'stacked','EdgeColor',colg(1,:),'LineWidth',1.5);
for i = 1:3
    set(h(i),'FaceColor',colg(i,:))
end
clear i h
h = barh(5,NA_smpbar_asv_frac(5,:),'stacked','EdgeColor',cols(1,:),'LineWidth',1.5);
for i = 1:3
    set(h(i),'FaceColor',cols(i,:))
end
clear i h
h = barh(6,NA_smpbar_asv_frac(6,:),'stacked','EdgeColor',colo(1,:),'LineWidth',1.5);
for i = 1:3
    set(h(i),'FaceColor',colo(i,:))
end
clear i h
axis ij
yticks([1:1:6])
yticklabels({'Diatoms','Dinoflagellates','Hacrobia','Chlorophytes','Dictyos+Pelagos','Other Ochrophyta'})
set(gca,'fontsize',20)
xlabel('Proportion of exported relative sequence abundance in surface samples')
axis([0 1 0 7])
title('North Atlantic')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Figure 5
%Relative sequence abundance of 6 major groups in each sample type by basin
%Calculate the relative abundance of 6 phytoplankton groups for each type
colors = [051 117 056;093 168 153;148 203 236;194 106 119;159 074 150;126 041 084]./255;

figure(5),clf
bbar = bar(all_frac,'stacked');
for i = 1:6
    set(bbar(i),'FaceColor',colors(i,:))
end
hold on
ylabel('Mean relative sequence abundance','fontsize',30)
xticklabels({'Bulk particles (35)','All other ind. particles (293)','Salp fecal pellets (27)','','Bulk particles (17)','Short pellets (56)','Aggregates + dense detritus (185)','Large loose + long fecal pellets (203)'})
xtickangle(35)
set(gca,'fontsize',24)
axis([0 9 0 1])
box on
legend({'Diatoms','Dinoflagellates','Hacrobia','Chlorophytes','Dictyos+Pelagos','Other Ochrophyta'},'fontsize',28)
clear i

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Figure 6
%EOFs with bulk ASVs shared with surface + POC flux
%Modeled POC flux vs. diatoms+hacrobia/all other phytos

%Isolate ASVs that are found in bulk and surface samples
%North Pacific
NP_taxall = NP_taxoP;

for i = 1:1519
    if NP_taxall.Genus(i)=='NaN' && NP_taxall.Class(i)=='NaN'
        allt_NP(i,1) = 'Other';
    elseif NP_taxall.Genus(i)=='NaN'
        allt_NP(i,1) = NP_taxall.Class(i);
    else allt_NP(i,1) = NP_taxall.Genus(i);
    end
end
clear i

NP_taxall.Class(NP_taxall.Class=='NaN') = 'Other';
NP_taxall.Class(NP_taxall.Division=='Dinoflagellata'&NP_taxall.Class=='Other') = 'Dinophyceae';
NP_taxall.Class(NP_taxall.Division=='Chlorophyta'&NP_taxall.Class=='Other') = 'Chlorophyta';
NP_taxall.Class(NP_taxall.Division=='Haptophyta'&NP_taxall.Class=='Other') = 'Haptophyta';

[~,indT] = sort(NP_taxall.Class);
NP_taxall = NP_taxall(indT,:);
allt_NP = allt_NP(indT,:);
NP_featall = NP_featP(indT,:);
clear indT

[a,~] = ismember(NP_taxall.FeatureID,bulktax_NP.FeatureID);
[b,~] = ismember(NP_taxall.FeatureID,parttax_NP.FeatureID);
[c,~] = ismember(NP_taxall.FeatureID,surftax_NP.FeatureID);

NPall_cat = categorical(1519,1);
for i = 1:1519
    if a(i)==0
        NPall_cat(i,1) = 'None';
    else
        NPall_cat(i,1) = NP_taxall.Genus(i);
    end
    if b(i)==0
        NPall_cat(i,2) = 'None';
    else
        NPall_cat(i,2) = NP_taxall.Genus(i);
    end
    if c(i)==0
        NPall_cat(i,3) = 'None';
    else
        NPall_cat(i,3) = NP_taxall.Genus(i);
    end
end
clear a b c

%NP: bulk, particles, surf wc
pall = strcmp(string(NPall_cat(:,[1,2,3])),'None');

%Repeat for North Atlantic

%6 major phytoplankton groups for both basins:
%NA_bulk_surf = separate all ASVs into 6 major groups
%NP_bulk_surf = separate all ASVs into 6 major groups

all_bulk_surf = [NP_bulk_surf;NA_bulk_surf];

%Load POC flux data
cd ~/EXPORTS/POC
load poc_std.mat %all_pocMs

%Calculate diatom+hacrobia/all other phytoplankton
diahac_all = (all_bulk_surf(:,1)+all_bulk_surf(:,3))./sum(all_bulk_surf(:,[2,4:6]),2);

%EOFs
EOF = [all_bulk_surf,all_pocMs];
EOFlabel = {'Diatoms','Dinoflagellates','Hacrobia','Chlorophytes','Dictyo+Pelago','Other Ochrophyta','POC flux'};

eofmeans = nanmean(EOF);
eofstd = nanstd(EOF);
for i = 1:size(EOF,1);
    eofs_center(i,:) = (EOF(i,:) - eofmeans)./eofstd;
end
clear i

[EOFs_all,AFs_all,eigvalues_all] = pca(eofs_center,'Centered',false,'Rows','complete'); 

var_explained = (eigvalues_all(1:end)/sum(eigvalues_all))';

%Plot figure 5A + B
figure(5),clf
subplot(1,2,1)
h2 = biplot(EOFs_all,'Scores',AFs_all(:,1:2),'VarLabels',plotlabel);
hold on
for k = 1:7
    h2(k).Color = 'k';
    h2(k).LineWidth = 2;
end
clear k
for k = 15:21
    h2(k).FontSize = 22;
end
clear k
for k = 22:50
    h2(k).MarkerSize = 55;
    h2(k).MarkerEdgeColor = [0.2734    0.5078    0.7031];
end
clear k
for k = 51:65
    h2(k).MarkerSize = 55;
    h2(k).MarkerEdgeColor = [0.755555555555556,0.226666666666667,0.528888888888889];
end
clear k
set(gca,'fontsize',30)
box on
legend([h2(22) h2(51)],{'North Pacific','North Atlantic'})
axis square

subplot(1,2,2)
hold on
errorbar(diahac_all(1:28),all_pocMs(1:28),poc_std(1:28),'LineStyle','none','Color',[0.2734    0.5078    0.7031],'LineWidth',2)
b1 = scatter(diahac_all(1:28),all_pocMs(1:28),350,repmat([0.2734    0.5078    0.7031],28,1),'filled');
errorbar(diahac_all(29:end),all_pocMs(29:end),poc_std(29:end),'LineStyle','none','Color',[0.755555555555556,0.226666666666667,0.528888888888889],'LineWidth',1.5)
b2 = scatter(diahac_all(29:end),all_pocMs(29:end),350,repmat([0.755555555555556,0.226666666666667,0.528888888888889],15,1),'filled');
b3 = plot(x,x*3.364+0.4709,'-k','LineWidth',4);
b4 = plot(x,x*1.17+0.97,':','Color',[0.2734    0.5078    0.7031],'LineWidth',3);
b5 = plot(x,x*3.4+0.24,'-.','Color',[0.755555555555556,0.226666666666667,0.528888888888889],'LineWidth',3);
l = legend([b1 b2 b3 b4 b5],{'North Pacific data','North Atlantic data','Linear fit (all data)','Linear fit (NP)','Linear fit (NA)'});
set(l,'fontsize',21)
xticks([0,.5,1,1.5,2,2.5,3])
xticklabels({'0.0','0.5','1.0','1.5','2.0','2.5','3.0'})
ylabel('Measured POC flux (mmol m^{-2} day^{-1})')
xlabel('(Diatoms+Hacrobia)/all other phytoplankton')
text(0.5,13,'y = 3.4x+0.47','fontsize',30)
text(0.5,12,'{\itR^2} = 0.72','fontsize',30)
text(0.5,11,'{\itP} value <<0.001','fontsize',30)
set(gca,'fontsize',30)
axis square
axis([0 3 0 15])
box on
clear b1 b2 b3 b4 b5


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Figure 7
%For each deployment, plot surface phytoplankton groups, bulk trap phytoplankton groups, and DH/all ratio
%Example shown with NA deployment 3: repeat for all
%Use depth metadata for each sample
%Use data for 6 surface groups and 6 bulk trap groups

figure(7)
subplot(3,4,11)
hold on
bbar = bar(35,surf6(6,:),70,'stacked');
hold on
for i = 1:6
    set(bbar(i),'FaceColor',colors(i,:))
end
clear i
a = area(NAdepth(12:16),deep6(12:16,:));
for i = 1:6
    a(i).FaceColor = colors(i,:);
    a(i).FaceAlpha = 0.9;
end
clear i
view(90,90)
text(400,0.6,'K','Color','w','fontsize',35)
xticks([0 100 200 300 400 500])
ylabel('Relative sequence abundance')
xlabel('Depth (m)')
set(gca,'fontsize',20)
axis([0 520 0 1])
box on

suplot(3,4,12)
t = tiledlayout(1,1);
ax1 = axes(t);
errorbar(all_pocMs(39:43),all_bulk_meta.Depth(39:43),poc_std(39:43),'horizontal','k','LineWidth',2)
hold on
scatter(ax1,all_pocMs(39:43),all_bulk_meta.Depth(39:43),300,'filled','k')
plot(ax1,all_pocMs(39:43),all_bulk_meta.Depth(39:43),'k','LineWidth',5)
axis ij
ax1.XColor = 'k';
set(gca,'fontsize',26,'box','off')
xlabel(ax1,'POC flux (mmol m^{-2} day^{-1})')
ylabel(ax1,'Depth (m)')
axis(ax1,[0 15 0 520])

ax2 = axes(t);
scatter(ax2,diahac_all(39:43),all_bulk_meta.Depth(39:43),200,repmat([0.755555555555556,0.226666666666667,0.528888888888889],5,1),'filled')
hold on
plot(ax2,diahac_all(39:43),all_bulk_meta.Depth(39:43),':','Color',[0.755555555555556,0.226666666666667,0.528888888888889],'LineWidth',5)
ax2.XColor = [0.755555555555556,0.226666666666667,0.528888888888889];
ax2.XAxisLocation = 'top';
ax2.Box = 'off';
axis ij
xlabel(ax2,'(Diatoms+Hacrobia)/all other phytoplankton')
set(gca,'fontsize',26,'Color','none')
axis(ax2,[0 2.75 0 520])
xline(ax2.XLim(2),'-k', 'linewidth',ax2.LineWidth);  
clear t ax1 ax2