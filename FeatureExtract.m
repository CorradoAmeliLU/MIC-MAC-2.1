%% Parameters

folder_path = '.\segmented\'; % Input folder (Imaris files).
folder_output = '.\extracted\'; % Output folder.

reprocess_samples = false; % If false, the script will not re-compute samples that have already been processed.

min_struct = 0; % Parameters to filter microglia by size (minimum 
max_struct = 10^8; % and maximum).

fill_width = 5;
se = strel('diamond', fill_width);

%% Opening Director

addpath(genpath('.\Library'));

dispstat('           -- MIC-MAC 2.0 -- ','keepthis');
dispstat('          -Feature Extraction- ','keepthis');
dispstat('          ','keepthis');

fold = dir(folder_path);
extr_fold = dir(folder_output);
[n_objs, ~] = size(fold);
samples = [""];
n_samps = 0;

for ii = 1 : n_objs
    if contains(fold(ii).name,'.ims') && ~contains(fold(ii).name,'.mat')
        if reprocess_samples
            n_samps = n_samps+1;
            samples(n_samps) = string(fold(ii).name);
        else
            n_occurrences = 0;
            for jj = 1:n_objs
                if(contains(fold(jj).name,fold(ii).name))
                    n_occurrences = n_occurrences + 1;
                end
            end
            if n_occurrences < 2
                n_samps = n_samps+1;
                samples(n_samps) = string(fold(ii).name);
            end
        end
    end
end

dispstat(['Found ' num2str(n_samps) ' samples.'],'timestamp','keepthis');

%% Extracting Sample by Sample

dispstat('Start Process. ','timestamp','keepthis');

metadata = table('Size',[n_samps,7],'VariableTypes',["string","double","double","double","double","double","double"],'VariableNames',["Sample","X","Y","Z","SUM","DENSITY","DENSITY_FILT"]);

for ii = 1: n_samps

    % Read Image
    FileIms = char(strcat(folder_path,samples(ii)))

    Image=bfopen(FileIms);
    Image = Image{1,1};

    [x,y]=size(Image{1,1});
    nStacks=length(Image);
    I=zeros(y,x,nStacks,'logical');

    for j=1:nStacks
        I(:,:,j)= Image{j,1}>0;
    end

    %Filling metadata information

    metadata.Sample(ii) = samples(ii);
    metadata.X(ii) = y;
    metadata.Y(ii) = x;
    metadata.Z(ii) = nStacks;
    metadata.SUM(ii) = sum(I,"all");
    metadata.DENSITY(ii) = metadata.SUM(ii)/(x*y*nStacks);

    %Filter structures by volume

    CC2 = bwlabeln(I);
    ST = regionprops3(CC2,'Volume');

    idx_between = find([ST.Volume] <= max_struct & [ST.Volume] > min_struct);
    I = ismember(CC2,idx_between);

    metadata.DENSITY_FILT(ii) = (sum(I,"all"))/(x*y*nStacks);

    % Retreiving Single Structures
    CC2 = bwlabeln(I);
    structs = regionprops3(CC2,'BoundingBox','Image');
    n_structs = height(structs);

    attributes = zeros(n_structs,43);

    for jj = 1 : n_structs

        dispstat(['         Computing Sample ' num2str(ii) ' out of ' num2str(n_samps) ': ' num2str(ceil(((jj-1)/n_structs)*100)) '%']);

        try

            % Isolating Structure
            bb = ceil(structs.BoundingBox(jj,:));
            vol = ismember(CC2, jj);
            vol = vol(bb(2):(bb(2)+bb(5)-1),bb(1):(bb(1)+bb(4)-1),bb(3):(bb(3)+(bb(6)-1)));

            [~,~,z] = size(vol);
            for kk = 1 : z
                vol(:,:,kk) = imclose(vol(:,:,kk),se);
            end

            % Computing Features
            [nodes_str,edges_str,~,~] = graph_formation(vol);
            [nodes_meas_table,graph_meas,~,~] = extract_graph_meas(nodes_str,edges_str);
            nodes_meas=table2array(nodes_meas_table);
            [id_node_center,~]=obtain_central_node(nodes_str,nodes_meas);
            [geom_meas,~]=extract_geom_meas(nodes_str,edges_str,vol,id_node_center);

            n_nodes = length(nodes_str);
            n_edges = length(edges_str);

            sum_edges = 0;
            for zz = 1 : n_edges
                if ~isempty(edges_str(zz).Length)
                    sum_edges = sum_edges + edges_str(zz).Length;
                end
            end

            attributes(jj,1:8) = cell2mat(struct2cell(geom_meas));
            attributes(jj,9:41) = cell2mat(struct2cell(graph_meas));

            attributes(jj,48) = n_nodes/attributes(jj,1);
            attributes(jj,49) = n_edges/attributes(jj,1);
            attributes(jj,50) = sum_edges/attributes(jj,1);

            %Rotating Cell to Compute Other Metrics

            max_ratio = 0;
            for angle = 1 : 5 : 90
                rotated_vol = imrotate3(vol,angle,[0 0 1],'nearest','loose','FillValues',0);
                [x,y,z] = size(rotated_vol);
                ratio = x/y;
                inv_ratio = y/x;
                if ratio > max_ratio
                    max_ratio = ratio;
                    correct_x = x;
                    correct_y = y;
                    correct_z = z;
                elseif inv_ratio > max_ratio
                    max_ratio = inv_ratio;
                    correct_x = y;
                    correct_y = x;
                    correct_z = z;
                end
            end

            attributes(jj,42) = correct_x;
            attributes(jj,43) = correct_y;
            attributes(jj,44) = correct_z;
            attributes(jj,45) = correct_x/correct_y;
            attributes(jj,46) = correct_x/(correct_y + (correct_z * 4.8));
            attributes(jj,47) = correct_y/(correct_z * 4.8);

        catch ME
        end

    end

    dispstat(['Sample ' samples(ii) ' Done.'],'timestamp','keepthis');

    % Export Data
    save(strcat(folder_output,samples(ii),'_attributes.mat'),'attributes')
    save(strcat(folder_output,samples(ii),'_volumes.mat'),'structs')
end

dispstat('End of Process. ','timestamp','keepthis');

writetable(metadata,'metadata2.csv')