
function transition_matrix = data_preprocess(file_path)
    graph_raw = importdata(file_path);
    data = graph_raw.data();
    %%% starting from 1 not 0, so adding 1 to sources and targets
    source_nodes = data(:,1)+1;
    target_nodes = data(:,2)+1;
    G = digraph(source_nodes,target_nodes);
    %%assign weights to edges with outgoing degrees
    for i = 1:length(source_nodes)
        out_degree = outdegree(G,source_nodes(i));
        weight(i) = 1/out_degree; 
    end
    transition_matrix = adjacency(G,weight);
end





