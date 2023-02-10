go_enrichment_analysis <- 
function (cluster_file = "clusters.txt", 
	  bk_gene_file = "back_ground_genes.txt", 
	  go_annotation_file = "data/Mouse.GO.annotation.txt", 
	  go_name_file = "data/GO.names.txt", 
	  gene_symbl_file = "data/Mouse.gene.symbl.txt")
{
	# Process GO annotation file. Only the GO terms with < 2500 gene annotations will be used for analysis.
	allgo = unique(read.table( go_annotation_file ))
	bk_genes = unique(readLines( bk_gene_file))
	total_gene_num = length(bk_genes)
	allgo = allgo[ allgo[,1] %in% bk_genes, ]
	bk_go_count = table(allgo[,2])
	idx = bk_go_count < 2500
	bk_go_count = bk_go_count[idx]
	allgo = allgo[allgo[,2] %in% attr(bk_go_count,"names"),]

	# Read in gene clusters
	clusters = unique( read.table( cluster_file ))
	cluster_size = table( clusters[,2] )
	cluster_size = sort( cluster_size, decreasing = T)
	cluster_name = attr( cluster_size, "names" )

	# Read in GO names
	go_name = read.table(go_name_file, sep="\t", comment.char = "", quote="")
	row.names(go_name) = go_name[,1]

	# Read in gene symbls
	gene_table = read.table(gene_symbl_file,sep="\t",comment.char = "", quote = "")
	gene_symbl = gene_table[,2]
	attr(gene_symbl,"names") = gene_table[,1]
	all_genes = unique(c(as.character(clusters[,1]),bk_genes))
	idx = all_genes %in% gene_table[,1]
	if ( sum(!idx) > 0){
		gene_symbl2 = all_genes[!idx]
		attr(gene_symbl2, "names") = gene_symbl2
		gene_symbl = c(gene_symbl, gene_symbl2)
	}

	# GO enrichment analysis
	for ( i in cluster_name){
		selected_genes = clusters[clusters[,2] == i, 1]
		cluster_size = length(selected_genes)
		cluster_go = allgo[allgo[,1] %in% selected_genes,]
		cluster_gene_symbl = gene_symbl[cluster_go[,1]]
		cluster_go_count = table(cluster_go[,2])
		go_ids = attr(cluster_go_count, "names")
		if ( length(go_ids) > 0){
			cluster_table = data.frame(cluster_id=character(), 
						   cluster_size=numeric(), 
						   go_rank = numeric(),
						   go_id=character(), 
						   go_catergory= character(), 
						   go_term=character(), 
						   cluster_go_count=numeric(),
						   genome_go_count=numeric(),
						   total_gene_number=numeric(),
						   genes_with_go_in_cluster=character(),
						   pValue=numeric())
			tt = 0
			for( j in go_ids){
				in_c = cluster_go_count[j]
				in_g = bk_go_count[j]
				pValue = phyper( in_c - 1, in_g, total_gene_num - in_g, cluster_size, lower.tail = F)
				idx = cluster_go[,2] == j
				genes_with_go = paste(sort(cluster_gene_symbl[idx]),collapse="/")
				tt = tt + 1
				cluster_table[tt,] <- c(i,cluster_size,1,j,go_name[j,2],go_name[j,3],in_c,in_g,total_gene_num,genes_with_go,pValue)
			}

			pValueAdjusted = p.adjust(cluster_table$pValue,"BH")
			cluster_table = cbind(cluster_table, pValueAdjusted)
			cluster_table = cluster_table[order(as.numeric(cluster_table$pValue)),]
			cluster_table$go_rank = 1:dim(cluster_table)[1]

			if (exists("all_table")){
				all_table = rbind(all_table, cluster_table)
			}else{
				all_table = cluster_table
			}
		}
	}

	if (exists("all_table")){
		idx = all_table$pValueAdjusted <= 0.05
		if( sum(idx) > 0){
			all_table = all_table[idx,]
			row.names(all_table) = NULL
			return( all_table )
		}else{
			cat( "No significant GO term found.","\n" )
			return( NULL )
		}
	}else{
		cat( "No significant GO term found.","\n" )
		return( NULL )
	}
}

mp_enrichment_analysis <- 
function (cluster_file = "clusters.txt", 
	  bk_gene_file = "back_ground_genes.txt", 
	  mp_annotation_file = "data/Mouse.MP.annotation.txt", 
	  mp_name_file = "data/MP.names.txt", 
	  gene_symbl_file = "data/Mouse.gene.symbl.txt")
{
	# Process MP annotation file. Only the MP terms with < 2500 gene annotations will be used for analysis.
	allmp = unique(read.table( mp_annotation_file ))
	bk_genes = unique(readLines( bk_gene_file))
	total_gene_num = length(bk_genes)
	allmp = allmp[ allmp[,1] %in% bk_genes, ]
	bk_mp_count = table(allmp[,2])
	idx = bk_mp_count < 2500
	bk_mp_count = bk_mp_count[idx]
	allmp = allmp[allmp[,2] %in% attr(bk_mp_count,"names"),]

	# Read in gene clusters
	clusters = unique( read.table( cluster_file ))
	cluster_size = table( clusters[,2] )
	cluster_size = sort( cluster_size, decreasing = T)
	cluster_name = attr( cluster_size, "names" )

	# Read in MP names
	mp_name = read.table(mp_name_file, sep="\t", comment.char = "", quote="")
	row.names(mp_name) = mp_name[,1]

	# Read in gene symbls
	gene_table = read.table(gene_symbl_file,sep="\t",comment.char = "", quote = "")
	gene_symbl = gene_table[,2]
	attr(gene_symbl,"names") = gene_table[,1]
	all_genes = unique(c(as.character(clusters[,1]),bk_genes))
	idx = all_genes %in% gene_table[,1]
	if ( sum(!idx) > 0){
		gene_symbl2 = all_genes[!idx]
		attr(gene_symbl2, "names") = gene_symbl2
		gene_symbl = c(gene_symbl, gene_symbl2)
	}

	# MP enrichment analysis
	for ( i in cluster_name){
		selected_genes = clusters[clusters[,2] == i, 1]
		cluster_size = length(selected_genes)
		cluster_mp = allmp[allmp[,1] %in% selected_genes,]
		cluster_gene_symbl = gene_symbl[cluster_mp[,1]]
		cluster_mp_count = table(cluster_mp[,2])
		mp_ids = attr(cluster_mp_count, "names")
		if ( length(mp_ids) > 0){
			cluster_table = data.frame(cluster_id=character(), 
						   cluster_size=numeric(), 
						   mp_rank = numeric(),
						   mp_id=character(),  
						   mp_term=character(),
						   mp_description=character(),
						   cluster_go_count=numeric(),
						   genome_go_count=numeric(),
						   total_gene_number=numeric(),
						   genes_with_mp_in_cluster=character(),
						   pValue=numeric())
			tt = 0
			for( j in mp_ids){
				in_c = cluster_mp_count[j]
				in_g = bk_mp_count[j]
				pValue = phyper( in_c - 1, in_g, total_gene_num - in_g, cluster_size, lower.tail = F)
				idx = cluster_mp[,2] == j
				genes_with_mp = paste(sort(cluster_gene_symbl[idx]),collapse="/")
				tt = tt + 1
				cluster_table[tt,] <- c(i,cluster_size,1,j,mp_name[j,2],mp_name[j,3],in_c,in_g,total_gene_num,genes_with_mp,pValue)
			}

			pValueAdjusted = p.adjust(cluster_table$pValue,"BH")
			cluster_table = cbind(cluster_table, pValueAdjusted)
			cluster_table = cluster_table[order(as.numeric(cluster_table$pValue)),]
			cluster_table$mp_rank = 1:dim(cluster_table)[1]

			if (exists("all_table")){
				all_table = rbind(all_table, cluster_table)
			}else{
				all_table = cluster_table
			}
		}
	}

	if (exists("all_table")){
		idx = all_table$pValueAdjusted <= 0.05
		if( sum(idx) > 0){
			all_table = all_table[idx,]
			row.names(all_table) = NULL
			return( all_table )
		}else{
			cat( "No significant MP term found.","\n" )
			return( NULL )
		}
	}else{
		cat( "No significant MP term found.","\n" )
		return( NULL )
	}
}

