package project.files;

import project.algorithm.HS.*;
import project.algorithm.NSGAII.*;
import project.data.*;

import java.io.*;
import java.util.*;

public class SaveFile {

	public SaveFile() {
		
	}
	
	public void saveHarmonyResults(String folder_name, int[] panel_index, ArrayList<Harmony> best_harmony) 
	{
		// Saving cut offs and pair (Sensitivity, Specificity) in a file
		String pre_name = folder_name;
		
		File file = new File(folder_name);
		if (!file.exists()) {
			file.mkdirs();
		}
		
		for (int i = 0; i < panel_index.length; i++) 
		{
			pre_name = pre_name+panel_index[i]+"-";
		}
		String cutoff_file = pre_name + "cutoff.csv";
		String var_file = pre_name + "value.csv";
		
		try
		{
			BufferedWriter bw_cut = new BufferedWriter(new FileWriter(cutoff_file));
			BufferedWriter bw_var = new BufferedWriter(new FileWriter(var_file));
			
			String content = "Sensitivity, Specificity";
			bw_var.write(content);
			bw_var.newLine();
			
			for (int i = 0; i < best_harmony.size(); i++) {
				bw_cut.write(Arrays.toString(best_harmony.get(i).getHarmonyValue()));
				bw_var.write(best_harmony.get(i).getSensitivity() + "," + best_harmony.get(i).getSpecificy());
				
				bw_cut.newLine();
				bw_var.newLine();
			}
			
			bw_cut.close();
			bw_var.close();
		} 
		catch(IOException ioe) 
		{
			ioe.printStackTrace();
		}
	}
	
	public void saveNSGAIIResults(String folder_name, String iterations, ArrayList<Individual> population, double time, ArrayList<Peptide> Peptides) 
	{
		String pre_name = folder_name+iterations+".csv";
		
		File file = new File(folder_name);
		if (!file.exists()) {
			file.mkdirs();
		}
		
		try 
		{
			BufferedWriter bw_value = new BufferedWriter(new FileWriter(pre_name));
			
			String content = "The running time is, " + time;
			bw_value.write(content);
			bw_value.newLine();
			
			int panel_size = population.get(1).getChromosomeSize();
			
			content = "";
			
//			for (int p = 1; p <= panel_size; p++)
//			{
//				content = content + "Peptide " + p + ',';
//			}
			
			content = content + "Panel, Sensitivity, Specificity, CutOffs, Rank, Fitness";
			bw_value.write(content);
			bw_value.newLine();
			
			// Saving results (Panel, Sensitivity, Specificity, Cutoffs, Rank and Fitness) per individual in the last population
			for (int i = 0; i < population.size(); i++) 
			{
				
				Individual chr = population.get(i);
				
				ArrayList<Integer> index = chr.getChomosomeIndex();
				ArrayList<Integer> index_copy = new ArrayList<Integer>();
				
				String peptide_names = "";
				for (int j = 0; j < index.size(); j++) 
				{
					int a = index.get(j)+1;
					//if (j == 0)
					  // peptide_names = peptide_names + Peptides.get(a).getPeptideName(); 
					//else
						//peptide_names = peptide_names + "," + Peptides.get(a).getPeptideName(); 
					index_copy.add(a);
					//System.out.println("\n" + peptide_names);
				}
				content = Arrays.toString(index_copy.toArray());
				content = content.replace(",","|");
				content = content + "," + chr.getSensitivity()+","+chr.getSpecificity()+",";
				double[] value = chr.getCutOff();
				double[] cut_off = new double[index.size()];
				for (int j = 0; j < cut_off.length; j++) {
					cut_off[j] = value[index.get(j)];
				}
				String cut = Arrays.toString(cut_off);
				cut = cut.replace(",","|");
				content = content+cut+ "," + chr.getRank() + "," + chr.getFitnessValue();
				
				bw_value.write(content);
				bw_value.newLine();
			}
			bw_value.close();
		} catch(IOException ioe) {
			ioe.printStackTrace();
		}
	}
	
	public void saveNSGAIIResults(int front_size, double time, String folder_name, String iterations, ArrayList<Individual> population, ArrayList<Peptide> Peptide_info, Global_GA_HS_Configurations set) 
	{
		String pre_name = folder_name+iterations+"_"+set.getPanelSize()+"_"+set.getThreshold()+".csv";
		
		File file = new File(folder_name);
		if (!file.exists()) {
			file.mkdirs();
		}
		
		try {
			BufferedWriter bw_value = new BufferedWriter(new FileWriter(pre_name));
			
			String content = "Experiment ID, Experiment Description, ";
			for (int i = 0; i < set.getPanelSize(); i++) {
				int index = i+1;
				content = content+"Peptide "+index+", ";
			}
			
			for (int i = 0; i < set.getPanelSize(); i++) {
				int index = i+1;
				content = content+"Cutoff "+index+", ";
			}
			content = content+"Sensitivity, Specificity, Panel Size, Threshold, Dataset Name";
			bw_value.write(content);
			bw_value.newLine();
			
			System.out.println(front_size);
			
			ArrayList<Individual> saved_chr = new ArrayList<Individual>();
			for (int i = 0; i < front_size; i++) {
				
				if (i > population.size()-1) {
					break;
				}
				
				Individual chr = new Individual(population.get(i).getChomosome().length);
				chr = population.get(i);
				if (SavedIndividual(chr, saved_chr)) {
					continue;
				} else {
					
					content = set.getExperimentDescription()+",";
					ArrayList<Integer> index = chr.getChomosomeIndex();
					for (int j = 0; j < index.size(); j++) {
						content = content+Peptide_info.get(index.get(j)).getPeptideName()+",";
					}
					content = fillEmpty(index.size(), set.getPanelSize(), content);
					
					double[] cut = chr.getCutOff();
					for (int j = 0; j < index.size(); j++) {
						content = content+cut[index.get(j)]+",";
					}
					content = fillEmpty(index.size(), set.getPanelSize(), content);
					
					content = content + chr.getSensitivity()+","+chr.getSpecificity()+","+set.getPanelSize()+","+set.getThreshold();
					
					bw_value.write(content);
					bw_value.newLine();
					
					saved_chr.add(chr);
				}
			}
			bw_value.close();
		} catch(IOException ioe) {
			ioe.printStackTrace();
		}
	}
	
	public boolean SavedIndividual(Individual chr, ArrayList<Individual> saved) 
	{
	
			for (int i = 0; i < saved.size(); i++) {
				Individual temp = saved.get(i);
				if (Arrays.equals(chr.getChomosome(), temp.getChomosome()) && chr.getSensitivity() == temp.getSensitivity() && chr.getSpecificity() == temp.getSpecificity()) {
					return true;
				}
			}
		
		return false;
	}
	
	public String fillEmpty(int panel_size, int max_panel_size, String content) {
		
		if (panel_size < max_panel_size) {
			int diff = max_panel_size - panel_size;
			while (diff > 0) {
				content = content + ",";
				diff--;
			}
		}
		
		return content;
		
	}
	
	public void saveGAResults(String folder_name, ArrayList<double[]> best_values, double time) 
	{
		String pre_name = folder_name+"paretofront.csv";
		
		File file = new File(folder_name);
		if (!file.exists()) {
			file.mkdirs();
		}
		
		try {
			BufferedWriter bw_value = new BufferedWriter(new FileWriter(pre_name));
			
			String content = "The running time is, " + time;
			bw_value.write(content);
			bw_value.newLine();
			
			content = "Sensitivity, Specificity";
			bw_value.write(content);
			bw_value.newLine();
			
			for (int i = 0; i < best_values.size(); i++) {
				content = best_values.get(i)[0]+","+best_values.get(i)[1];
				bw_value.write(content);
				bw_value.newLine();
			}
			bw_value.close();
		} catch(IOException ioe) {
			ioe.printStackTrace();
		}
	}
}
