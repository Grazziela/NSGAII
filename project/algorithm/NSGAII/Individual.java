package project.algorithm.NSGAII;

import project.algorithm.HS.Harmony;
import project.data.*;
import project.start.*;
import project.files.*;

import java.io.*;
import java.util.*;

public class Individual 
{
	private int[] chromosome = null;
	private int totalgen = 0;
	private ArrayList<Integer> gen = null;
	private ArrayList<Integer> nogen = null;
	private double fitness_value;
	private double sum_sen_spe;

	private int rank = 0;
	private int number_dominates = 0;
	private ArrayList<Integer> dominators = null;
	private double crowd_distance = 0.0;

	private Harmony har = null;

	public Individual() 
	{
	}

	public Individual(int length) 
	{
		this.chromosome = new int[length];
		this.har = new Harmony(length);

		clear();
	}

	public void buildIndividual(int panel_size, RandomGenerator random_generator, ArrayList<Sample> sample_values) 
	{
		int[] chromosome = new int[this.chromosome.length];
		int flag = 0;
		int index = 0;
		
		// This is an array of indexes from 1..numberPeptides to avoid that the random search points to the same peptide twice in the panel
		ArrayList<Integer> nogen_copy = (ArrayList<Integer>)this.nogen.clone();
		
		

		for (int i = 0; i < panel_size; i++) 
		{
			//flag = 0;
			//while(flag == 0) 
			//{
				// Generating a random values from 0 to numberPeptides - 1
				index = (int) random_generator.getRandomIndex(nogen_copy.size());
								
				//if(sample_values.get(sample_values.size()-1).getPeptideValue(index) != 0.0)
			//	{
				//	flag = 1;
					chromosome[nogen_copy.get(index)] = 1;
					nogen_copy.remove(nogen_copy.get(index));
			//	}

			//}
		}
		buildIndividual(chromosome);
	}

	public void buildIndividual(int[] chromosome)
	{
		this.chromosome = new int[chromosome.length];
		clear();
		for (int i = 0; i < chromosome.length; i++) 
		{
			this.chromosome[i] = chromosome[i];
			if (chromosome[i] == 1) 
			{
				this.totalgen++;
				this.gen.add(i);
				this.nogen.remove((Integer) i);
			}
		}
	}

	/*
	 * flip the gen value from 1 to 0.
	 */
	public void flipGen(int size, int random_seed) 
	{
		RandomGenerator rnd = new RandomGenerator(random_seed);
		while(this.totalgen > size) {
			int index = (int)rnd.getRandomIndex(this.totalgen);
			int val = this.gen.get(index);
			this.chromosome[val] = 0;
			this.gen.remove(index);
			this.nogen.add(val);
			this.totalgen--;
		}
	}

	/*
	 * flip the gen value from 0 to 1.
	 */
	public void antiflipGen(int size, int random_seed) 
	{
		RandomGenerator rnd = new RandomGenerator(random_seed);
		while (this.totalgen < size) {
			int index = (int)rnd.getRandomIndex(this.nogen.size());
			int val = this.nogen.get(index);
			this.chromosome[val] = 1;
			this.nogen.remove(index);
			this.gen.add(val);
			this.totalgen++;
		}
	}
	
	public int getChromosomeSize()
	{
		return this.chromosome.length;
	}

	public double getFitnessValue() 
	{
		return this.fitness_value;
	}


	public double getSumValue() 
	{
		return this.sum_sen_spe;
	}

	public double getSensitivity() 
	{
		return this.har.getSensitivity();
	}

	public double getSpecificity() 
	{
		return this.har.getSpecificy();
	}

	public int getNumberDominates() 
	{
		return this.number_dominates;
	}

	public void updateNumberDominates() 
	{
		this.number_dominates++;
	}

	public void decreaseNumberDominates() 
	{
		this.number_dominates--;
	}

	public void setNumberDominates(int number_dominates) 
	{
		this.number_dominates = number_dominates;
	}

	public int getRank() 
	{
		return this.rank;
	}

	public void setRank(int rank) 
	{
		this.rank = rank;
	}

	public double getCrowdDistance() 
	{
		return this.crowd_distance;
	}

	public void setCrowdDistance(double crowd_distance) 
	{
		this.crowd_distance = crowd_distance;
	}

	public ArrayList<Integer> getDominators() 
	{
		return this.dominators;
	}

	public void setDominators(ArrayList<Integer> dominators) 
	{
		this.dominators.clear();
		for (int i = 0; i < dominators.size(); i++)
		{
			int value = dominators.get(i);
			this.dominators.add(value);
		}
	}

	public int getSpecificDominator(int index) 
	{
		return this.dominators.get(index);
	}

	public void removeSpecificDominator(int index) 
	{
		this.dominators.remove(index);
	}

	public double[] getCutOff() 
	{
		return this.har.getHarmonyValue();
	}

	public void setFitnessValue(double fitness_value) 
	{
		this.fitness_value = fitness_value;
	}

	public void setSumValues(double sum_sen_spe) 
	{
		this.sum_sen_spe = sum_sen_spe;
	}

	public void setFullChomosome(int[] chromosome) 
	{
		buildIndividual(chromosome);
	}

	public void setChomosomeValue(int index, int value) 
	{
		if (this.chromosome[index] != value)
		{
			if (value == 1)
			{
				this.totalgen++;
				this.gen.add(index);
				this.nogen.remove((Integer) index);
			} 
			else 
			{
				this.totalgen--;
				this.nogen.add(index);
				this.gen.remove((Integer) index);
			}
		}

		this.chromosome[index] = value;
	}

	public int getChomosomeValue(int index) 
	{
		return this.chromosome[index];
	}

	public int[] getChomosome() 
	{
		int[] cloned = new int[this.chromosome.length];
		for (int i = 0; i < cloned.length; i++) 
		{
			cloned[i] = this.chromosome[i];
		}

		return cloned;
	}

	public void setHarmony(Harmony harmony) 
	{
		this.har.setHarmony(harmony, this.chromosome);
		setFitnessValue(harmony.getFitness());
		double sensitivity = harmony.getSensitivity();
		double specificity = harmony.getSpecificy();
		double sum = sensitivity + specificity;
		setSumValues(sum);
	}

	public ArrayList<Integer> getChomosomeIndex() 
	{
		return (ArrayList<Integer>)this.gen.clone();
	}

	public int getTotalGen() 
	{
		return this.totalgen;
	}

	public void clear() 
	{
		this.chromosome = new int[this.chromosome.length];
		this.totalgen = 0;
		this.fitness_value = Double.MAX_VALUE;
		this.sum_sen_spe = Double.MAX_VALUE;
		this.gen = new ArrayList<Integer>();
		this.nogen = new ArrayList<Integer>();
		for (int i = 0; i < this.chromosome.length; i++) {
			this.nogen.add(i);
		}

		NSGAIIInitialisation();
	}

	public void NSGAIIInitialisation() 
	{
		this.rank = 0;
		this.dominators = new ArrayList<Integer>();
		this.number_dominates = 0;
		this.crowd_distance = 0.0;
	}

	public void PrintOnScreen() 
	{
		System.out.println("The fitness value is " + this.fitness_value + " and sensitivity is "+this.har.getSensitivity() + " and the specificity is "+ this.har.getSpecificy());
		System.out.println("The panel is " + Arrays.toString(this.gen.toArray()));
		System.out.println("The cut-off of this panel is "+ Arrays.toString(this.har.getHarmonyValue()));
	}

	public static Comparator<Individual> sortByFITNESS = new Comparator<Individual>() 
	{
		public int compare(Individual ch1, Individual ch2) 
		{

			double fitness1 = ch1.getFitnessValue();
			double fitness2 = ch2.getFitnessValue();

			/*For ascending order*/
			if (fitness1-fitness2 > 0) 
			{
				return 1;
			} 
			if (fitness1-fitness2 < 0)
			{
				return -1;
			}
			return 0;  
		}
	};

	public static Comparator<Individual> sortBySENS = new Comparator<Individual>() 
	{
		public int compare(Individual ch1, Individual ch2) {

			double sens1 = ch1.getSensitivity();
			double sens2 = ch2.getSensitivity();

			/*For descending order*/
			if (sens1-sens2 > 0) {
				return -1;
			} 

			if (sens1-sens2 < 0){
				return 1;
			} 

			return 0;
		}
	};   

	public static Comparator<Individual> sortBySPEC = new Comparator<Individual>() {

		public int compare(Individual ch1, Individual ch2) {

			double spec1 = ch1.getSpecificity();
			double spec2 = ch2.getSpecificity();

			/*For descending order*/
			if (spec1-spec2 > 0) {
				return -1;
			} 

			if (spec1-spec2 < 0){
				return 1;
			} 

			return 0;
		}
	}; 	   

	public static Comparator<Individual> sortByDISTANCE = new Comparator<Individual>() 
	{
		public int compare(Individual ch1, Individual ch2) {

			double distance1 = ch1.getCrowdDistance();
			double distance2 = ch2.getCrowdDistance();

			/*For descending order*/
			if (distance1-distance2 > 0) {
				return -1;
			} 

			if (distance1-distance2 < 0){
				return 1;
			} 

			return 0;
		}
	}; 				   

	@Override
	public Individual clone() 
	{
		Individual cloned = null;
		try 
		{
			cloned = (Individual) super.clone();
		} 
		catch (CloneNotSupportedException e) 
		{

		}

		return cloned;
	}
}
