package project.start;

import project.algorithm.*;
import project.data.*;
import project.files.*;

import java.util.*;
import java.io.*;

/*
 Main file. This is where the data processing and search for the solutions start.
 Run this module after altering the file name to the file you want to create the cut offs/establish the screening panels for.
 This code solve the problems of determining a pareto front for possible panels of Peptides/peptides/genes and their respective 
 cut offs that separate/categorise the dataset instances from tests and controls (wild types)
 The NSGA algorithm searches for the best panels of peptides/antigens etc in the Pareto front
 The Harmony Search establishes the cut off values for each peptide within a panel.
 
  Note that all configurations for the GA and HS are in the class file Global_GA_HS_Configurations. 
  So if parameters need to be adjusted, that file should be modified. 
  
  In the next development iteration, those parameters will be read from a file, to allow for better control 
  of multiple experiments and their documentation.
 */
public class ProjectNSGAII {

	public static void main(String[] args) {

		/* data file, where the information about samples (rows) and peptides (columns) will be available.
		 * The code accepts a .csv file format, where the last column is the class of each sample: test (value = 1) or control (value = 0)
		*/
		String input_data_file_name = ".\\Top280Peptides_TrainingSetSize_23_Samples_AllData.csv";
		String experiment_description = "Peptides Experiment - Training set size 23";

		// Reading input data
		ReadInputDatasetFile data = new ReadInputDatasetFile(input_data_file_name);

		int max_panel_size = 100; // maximum size of the panel within the Pareto Front that satisfies the constraints of max(Sensitivity, Specificity)	
		
		// Global_GA_HS_Configurations is a separate class containing all global settings for the GA, HS and files used. 
		Global_GA_HS_Configurations global_configurations = new Global_GA_HS_Configurations(1, data.getChromosomeSize(), experiment_description);


		// This loop will iterate over many panel sizes, from n0 to nk and increment j. So solutions for multiple panels will be generated and saved 
		for (int i = 20; i < max_panel_size; i = i+5) 
		{
			global_configurations.setPanelSize(i);
			
			// Finding the Pareto Front with resulting panels and saving the results
			String save_folder = ".\\PanelSize "+ i +"\\";
			String final_result_folder = ".\\results\\";
			
			Solver solver = new Solver(i);
			solver.start(data, save_folder, final_result_folder, global_configurations);
		}
	}

	
	// STILL NEED TO FIGURE OUT WHAT THE METHODS BELOW DO, BUT i SUSPECT THEY HAVE BEEN CREATED TO PUT THE FILE IN THE FORMAT FOR TABLEAU VIZ
	public static void testFunction() {

		String file_name = "C:\\Users\\uizch1\\OneDrive - The University of Nottingham\\NSGAIIProject-master\\src\\project\\TEST.csv";
		String save_folder = "C:\\Users\\uizch1\\OneDrive - The University of Nottingham\\NSGAIIProject-master\\src\\project";

		int[] panel_size_set = {6, 8, 10, 12, 15};
		int[] threshold_set = {1};

		File file = new File(save_folder);
		if (!file.exists()) {
			file.mkdirs();
		}

		ReadInputDatasetFile rf = new ReadInputDatasetFile(file_name);

		ArrayList<Peptide> Peptide_information = new ArrayList<Peptide>();
		Peptide_information = (ArrayList<Peptide>)rf.getPeptideNames().clone();

		file_name = "C:\\Users\\uizch1\\OneDrive - The University of Nottingham\\NSGAIIProject-master\\src\\project\\data.csv";

		try 
		{
			BufferedReader br = new BufferedReader(new FileReader(file_name));
			String information = "";
			ArrayList<String[]> information_set = new ArrayList<String[]>();

			ArrayList<ArrayList<ArrayList<String[]>>> information_collection = new ArrayList<ArrayList<ArrayList<String[]>>>();

			int panel_size = -1;
			int threshold = -1;
			String time = "";
			while ((information = br.readLine()) != null) {
				if (!information.trim().equals("")) {
					if (information.contains("*")) {
						// PrintOutInformation(information_set, panel_size, threshold, time, save_folder, Peptide_information);
						if (threshold - 1 >= information_collection.size()) {
							ArrayList<ArrayList<String[]>> temp = new ArrayList<ArrayList<String[]>>();
							temp.add((ArrayList<String[]>)information_set.clone()); 
							information_collection.add(temp);
						} 
						else 
						{
							ArrayList<ArrayList<String[]>> temp = new ArrayList<ArrayList<String[]>>();
							temp = information_collection.get(threshold-1);
							temp.add((ArrayList<String[]>)information_set.clone());
							information_collection.remove(threshold-1);
							information_collection.add(threshold-1, temp);
						}
						information_set.clear();
					} 
					else 
					{
						if (information.contains("CutOff")) 
						{
							continue;
						} 
						else if (information.contains("panel")) {
							String[] info_set = information.split(",");
							String info = info_set[3].trim();
							panel_size = Integer.parseInt(info);
							info = info_set[5].trim();
							threshold = Integer.parseInt(info);
							time = info_set[1].trim();
						} 
						else 
						{
							String[] info = information.split(",");
							String[] copy = new String[info.length];
							for (int i = 0; i < info.length; i++) {
								if (i < panel_size) {
									copy[i] = info[i];
								} else if (i < panel_size+2) {
									copy[i+panel_size] = info[i];
								} else {
									copy[i-2] = info[i];
								}
							}
							information_set.add(copy);
						}
					}
				}
			}

			PrintOutInformation(information_collection, panel_size_set, threshold_set, save_folder, Peptide_information);
		} 
		catch(IOException e) 
		{
			System.out.println(e.getMessage());
			System.exit(-1);
		}
	}

	public static void PrintOutInformation(ArrayList<ArrayList<ArrayList<String[]>>> data, int[] panel_size, int[] threshold_size, String save_folder, ArrayList<Peptide> Peptide_info) {
		String file_name = save_folder+"result summary.csv";
		String save_info = "";
		int line_number = 0;

		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(file_name, true));

			for (int i = 0; i < data.size(); i++) {
				ArrayList<ArrayList<String[]>> temp = new ArrayList<ArrayList<String[]>>();
				temp = data.get(i);
				save_info = ""+threshold_size[i]+"PeptidePOSITIVEFORCANCERPOSITIVE";
				bw.write(save_info);
				bw.newLine();

				for (int j = 0; j < temp.size(); j++) {
					int diff = panel_size.length - temp.size();
					ArrayList<String[]> sol_information = temp.get(j);
					save_info = ", ";
					for (int k = 0; k < panel_size[j+diff]; k++) {
						int index = k+1;
						save_info = save_info+"Peptide "+index+", ";
					}

					for (int k = 0; k < panel_size[j+diff]; k++) {
						int index = k+1;
						save_info = save_info+"cutt off "+index+",";
					}

					save_info = save_info+"sensitivity, specificity";
					line_number++;

					bw.write(save_info);
					bw.newLine();

					for (int k = 0; k < sol_information.size(); k++) {
						save_info = line_number+", ";
						String[] info = new String[sol_information.get(k).length];
						info = sol_information.get(k).clone();

						for (int m = 0; m < info.length; m++) {
							if (m < panel_size[j+diff]) {
								String in = info[m].trim();
								int index = Integer.parseInt(in);
								save_info = save_info + Peptide_info.get(index-1).getPeptideName()+", ";
							} else {
								// System.out.println("panel size is "+panel_size+" and the related threshold is "+threshold);
								String in = info[m].trim();
								save_info = save_info + in + ", ";
							}
						}

						bw.write(save_info);
						bw.newLine();
						line_number++;
					}
					line_number = 0;

					bw.write(" ");
					bw.newLine();

				}
				line_number = 0;
				bw.write(" ");
				bw.newLine();
			}

			bw.close();
		} catch(IOException e) {
			System.out.println(e.getMessage());
			System.exit(-1);
		}
	}

	public static void PrintOutInformation(ArrayList<String[]> data, int panel_size, int threshold, String time, String save_folder, ArrayList<Peptide> Peptide_info) {
		String file_name = save_folder+"NewData_"+panel_size+"_"+threshold+".csv";
		// String file_name = save_folder+"result summary.csv";
		String save_info = "";
		int line_number = 0;

		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(file_name));

			save_info = ", ";
			for (int i = 0; i < panel_size; i++) {
				int index = i+1;
				save_info = save_info+"Peptide "+index+", ";
			}

			for (int i = 0; i < panel_size; i++) {
				int index = i+1;
				save_info = save_info+"cutt off "+index+",";
			}

			save_info = save_info+"sensitivity, specificity";
			line_number++;

			bw.write(save_info);
			bw.newLine();

			for (int i = 0; i < data.size(); i++) {
				save_info = line_number+", ";
				String[] info = new String[data.get(i).length];
				info = data.get(i).clone();

				for (int j = 0; j < info.length; j++) {
					if (j < panel_size) {
						String in = info[j].trim();
						int index = Integer.parseInt(in);
						save_info = save_info + Peptide_info.get(index-1).getPeptideName()+", ";
					} else {
						// System.out.println("panel size is "+panel_size+" and the related threshold is "+threshold);
						String in = info[j].trim();
						save_info = save_info + in + ", ";
					}
				}

				bw.write(save_info);
				bw.newLine();
				line_number++;
			}

			bw.write(" ");
			bw.newLine();

			bw.close();
		} catch(IOException e) {
			System.out.println(e.getMessage());
			System.exit(-1);
		}
	}
}
