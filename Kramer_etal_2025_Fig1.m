%Sasha Kramer
%MBARI

%%%EXPORTS 18S Aitchison clr transform + PCA
cd ~/18S/QIIME/Jan23
load EXPORTS_18S_20230106.mat

%Isolate samples to use
feat = table2array(feature_table(:,[1:311,400:886,893:1066,1083:1190]));
meta = meta_table([1:311,400:886,893:1066,1083:1190],:);
clear feature_table meta_table

taxo = taxo_table;
taxo(sum(feat,2)==0,:) = [];
feat(sum(feat,2)==0,:) = [];
clear taxo_table

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

taxo_p = taxo(indP==1,:);
feat_p = feat(indP==1,:);
clear feat taxo photo hetero indP

%Need to filter the data to remove rows with too many zeros
%Remove ASVs that are <2% of total read counts across all samples
filtfeat = feat_p;

filtfeat = filtfeat./sum(filtfeat);
ind = find(max(filtfeat,[],2)<0.02);

filtfeat = feat_p;
filtfeat(ind,:) = [];
filtfeatN = filtfeat./sum(filtfeat);
taxo_filt = taxo_p;
taxo_filt(ind,:) = [];
clear ind

%Compiling data at a higher taxo level
[~,ind] = sortrows(taxo_filt,[2,3,4,5,6,7,8,9]);

taxo_filt = taxo_filt(ind,:);
filtfeat = filtfeat(ind,:);
filtfeatN = filtfeatN(ind,:);
clear ind

taxo_filtT = taxo_filt(:,2:8);

%Index by genus
[U,I] = unique(taxo_filtT,'first');
U(61:70,:) = [];
I(61:70,:) = [];

for j = 1:979
    for i = 1:122
        if taxo_filtT.Class(j)==U.Class(i)&&taxo_filtT.Order(j)==U.Order(i)&&taxo_filtT.Family(j)==U.Family(i)&&taxo_filtT.Genus(j)==U.Genus(i)
            feat_ind(j,:) = i;
        end
    end
end
clear i j

for i = 1:122
    if length(find(feat_ind==i))>1
        feat_sum(i,:) = sum(filtfeat(feat_ind==i,:));
    elseif length(find(feat_ind==i))==1
        feat_sum(i,:) = filtfeat(feat_ind==i,:);
    end
end
clear i

feat_sumN = feat_sum./sum(feat_sum);

%Have to add 1 for pseudocount
feat1 = filtfeat+1;
featN1 = filtfeatN+1;
feat_sum1 = feat_sum+1;
feat_sumN1 = feat_sumN+1;

feat1_clr = log(feat1./repmat(geomean(feat1,2),1,size(feat1,2)));
featN1_clr = log(featN1./repmat(geomean(featN1,2),1,size(featN1,2)));
featf1_clr = log(feat_sum1./repmat(geomean(feat_sum1,2),1,size(feat_sum1,2)));
featfN1_clr = log(feat_sumN1./repmat(geomean(feat_sumN1,2),1,size(feat_sumN1,2)));

%Sort samples by type
yind = [1:166,312:348,799:917,167:201,973:1026,349:798,918:934,1027:1080];

%Perform PCA analysis
[EOFs_featf1,AFs_featf1,eigvalues_featf1] = svd(featf1_clr,'Centered',false,'Rows','complete'); 

var_explained = (eigvalues_featf1(1:end)/sum(eigvalues_featf1))';

%Eigenvalues/loadings
EOFs_f1fp = EOFs_featf1(:,1:2);

%Eigenvectors
AFs_f1fp = AFs_featf1(:,1:2);

%Define colors
newcolors = [051 117 056;093 168 153;148 203 236;194 106 119;159 074 150;126 041 084]./255;

ac = hex2rgb('#e09cc3');
ac(2,:) = hex2rgb('#ac395b');
ac(3,:) = hex2rgb('#dcaae8');

pc = hex2rgb('#97b2c9');
pc(2,:) = hex2rgb('#50738a');
pc(3,:) = hex2rgb('#b7e2d6');

nmcol2(1:322,:) = repmat(pc(1,:),322,1);
nmcol2(323:357,:) = repmat(pc(2,:),35,1);
nmcol2(358:411,:) = repmat(pc(3,:),54,1);
nmcol2(412:861,:) = repmat(ac(1,:),450,1);
nmcol2(862:878,:) = repmat(ac(2,:),17,1);
nmcol2(879:932,:) = repmat(ac(3,:),54,1);

%Make figures
figure(1),clf
h2 = biplot(AFs_f1fp,'VarLabels','');
hold on
for k = 1:45
    h2(k).Color = newcolors(2,:);
    h2(k).LineWidth = 2.5;
end
clear k
for k = 46:60
    h2(k).Color = newcolors(4,:);
    h2(k).LineWidth = 2.5;
end
clear k
for k = 61:78
    h2(k).Color = newcolors(3,:);
    h2(k).LineWidth = 2.5;
end
clear k
for k = 79:100
    h2(k).Color = newcolors(1,:);
    h2(k).LineWidth = 2.5;
end
clear k
for k = [101:110,114:116]
    h2(k).Color = newcolors(6,:);
    h2(k).LineWidth = 2.5;
end
clear k
for k = [111:113,117:122]
    h2(k).Color = newcolors(5,:);
    h2(k).LineWidth = 2.5;
end
clear k
for k = 122+(1:45)
    h2(k).MarkerEdgeColor = newcolors(2,:);
    h2(k).MarkerSize = 45;
end
clear k

for k = 122+(46:60)
    h2(k).MarkerEdgeColor = newcolors(4,:);
    h2(k).MarkerSize = 45;
end
clear k
for k = 122+(61:78)
    h2(k).MarkerEdgeColor = newcolors(3,:);
    h2(k).MarkerSize = 45;
end
clear k
for k = 122+(79:100)
    h2(k).MarkerEdgeColor = newcolors(1,:);
    h2(k).MarkerSize = 45;
end
clear k
for k = 122+([101:110,114:116])
    h2(k).MarkerEdgeColor = newcolors(6,:);
    h2(k).MarkerSize = 45;
end
clear k
for k = 122+([111:113,117:122])
    h2(k).MarkerEdgeColor = newcolors(5,:);
    h2(k).MarkerSize = 45;
end
clear k
set(gca,'fontsize',18)
set(gca,'fontsize',30)
xticks([-80 -40 0 40 80])
yticks([-80 -40 0 40 80])
xlabel('PC1 (26.5%)')
ylabel('PC2 (11.3%)')
box on
l = legend([h2(201) h2(123) h2(183) h2(168) h2(233) h2(223)],{'Diatoms','Dinoflagellates','Hacrobia','Chlorophytes','Dictyos+Pelagos','Other Ochrophyta'});
set(l,'fontsize',24)
axis square
clear l

figure(2),clf
h2 = biplot(AFs_f1fp,'Scores',EOFs_f1fp(yind,:),'VarLabels','');
hold on
for k = 1:122
    h2(k).Color = [0.65 0.65 0.65];
end
clear k
for k = 123:245
    h2(k).MarkerSize = 12;
    h2(k).MarkerEdgeColor = [0.65 0.65 0.65];
end
clear k
for k = [246:567,656:1106] 
    h2(k).Marker = 'o';
    h2(k).MarkerSize = 12;
    h2(k).LineWidth = 1;
    h2(k).MarkerEdgeColor = nmcol2(k-245,:);
end
clear k
for k = [568:655,1107:1177]
    h2(k).MarkerSize = 50;
    h2(k).MarkerEdgeColor = nmcol2(k-245,:);
end
clear k
set(gca,'fontsize',30)
xticks([-80 -40 0 40 80])
yticks([-80 -40 0 40 80])
xlabel('PC1 (26.5%)')
ylabel('PC2 (11.3%)')
box on
l = legend([h2(256) h2(579) h2(615) h2(669) h2(1119) h2(1136)],{'NP ind. particles','NP bulk particles','NP surface water','NA ind. particles','NA bulk particles','NA surface water'});
set(l,'fontsize',24)
axis square
clear l