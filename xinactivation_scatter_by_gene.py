'''
Purpose: scatter plot for escape, inactivated, and variable genes, 
    with genes on x-axis, average methylation on y-axis,
    points colored on whether they are XX placenta, XY placenta,
    decidua with XX placenta, decidua with XY placenta  
Input: tables of average methylation in gene body coordinates or promoter coordinates 
    for genes with x-inactivation status
Output: series of scatter plots for x-inactivation status gene groups
    
'''
import matplotlib.pyplot as plt

# input tables
xinact_gene_file = "/data/CEM/wilsonlab/projects/placenta/2021_placenta_test/methylation3/xinactivated_genes/XCI_data_all_samples_gene.txt"
xinact_gene_file2 = "/data/CEM/wilsonlab/projects/placenta/2021_placenta_test/methylation3/xinactivated_genes/XCI_data_placenta_samples_gene.txt"
xinact_gene_file3 = "/data/CEM/wilsonlab/projects/placenta/2021_placenta_test/methylation3/xinactivated_genes/XCI_data_decidua_samples_gene.txt"
xinact_promoter_file = "/data/CEM/wilsonlab/projects/placenta/2021_placenta_test/methylation3/xinactivated_genes/XCI_data_all_samples_promoter.txt"
xinact_promoter_file2 = "/data/CEM/wilsonlab/projects/placenta/2021_placenta_test/methylation3/xinactivated_genes/XCI_data_placenta_samples_promoter.txt"
xinact_promoter_file3 = "/data/CEM/wilsonlab/projects/placenta/2021_placenta_test/methylation3/xinactivated_genes/XCI_data_decidua_samples_promoter.txt"
data_files = [xinact_gene_file, xinact_gene_file2, xinact_gene_file3, xinact_promoter_file, xinact_promoter_file2, xinact_promoter_file3]

# hard-code needed info
placenta_XX_string = 'Placenta_XX'
placenta_XX_color = 'red'

placenta_XY_string = 'Placenta_XY'
placenta_XY_color = 'green'

decidua_w_XX_string = 'Decidua_XX'
decidua_w_XX_color = 'purple'

decidua_w_XY_string = 'Decidua_XY'
decidua_w_XY_color = 'gold'

xci_status_groups = ["escape","inactivated","variable"]

gene_index = 0
status_index = 1
data_start_index = 2

# iterate both tables
for file in data_files: 
    # load data
    inputhandle = open(file, 'r')
    data = inputhandle.readlines()
    inputhandle.close()

    headerline = data[0].replace('\n','')
    headers = headerline.split('\t')

    # iterate gene types
    for xci_status in xci_status_groups:

        # creating the points to plot
        x_coordinates = []
        y_coordinates = []
        gene_names = []
        point_colors = []

        # iterate through the lines of the data file creating points to plot
        counter = 0
        for line in data[1:len(data)]:
            line = line.replace("\n","")
            linesplit = line.split('\t')
            
            # go through data for each line collecting methylation for this XCI status
            for i in range(data_start_index,len(linesplit)):
                if (linesplit[status_index] == xci_status):
                    x_coordinates.append(counter)
                    y_coordinates.append(float(linesplit[i]))
                    if (not linesplit[gene_index] in gene_names):
                        gene_names.append(linesplit[gene_index])
                    if (placenta_XX_string in headers[i]):
                        point_colors.append(placenta_XX_color)
                    elif (placenta_XY_string in headers[i]):
                        point_colors.append(placenta_XY_color)
                    elif (decidua_w_XX_string in headers[i]):
                        point_colors.append(decidua_w_XX_color)
                    elif (decidua_w_XY_string in headers[i]):
                        point_colors.append(decidua_w_XY_color)
                    
            # after going through all the data for this gene, move on to next gene 
            if (len(x_coordinates) >= 1):
                counter = x_coordinates[-1] + 1

        #print(x_coordinates)
        #print(y_coordinates)
        #print(point_colors)
        #print(gene_names)

        # create plot
        w = 6
        if (xci_status != "escape"):
            w = 11

        fig,ax = plt.subplots(figsize = (w,4.8))

        index_list = [*range(0, len(gene_names), 1)] 

        ax.scatter(x = x_coordinates, y = y_coordinates, marker = "o", s = 14, c = point_colors, alpha=0.5)
        
        ax.set_ylabel("average percent methylation")
        ax.set_ylim(0,100)
        ax.set_xlabel("gene")
        title_string = xci_status + " genes ("
        if ("gene.txt" in file):
            title_string += "gene body coordinates)"
        elif ("promoter.txt" in file):
            title_string += "promoter coordinates)"
        ax.set_title (title_string)
        ax.set_xticks(index_list)
        ax.set_xticklabels(gene_names, rotation = 90, fontsize = 'x-small')

        # save figure
        
        fig.tight_layout(pad=2.0)
        outputpng = file.replace(".txt","_scatter_" + xci_status + ".pdf") 
        plt.savefig(outputpng)
        print ("Created plot: ", outputpng)

