package project.start;

import java.util.*;

public class RandomGenerator 
{
	private Random rnd = null;
	
	public RandomGenerator(int random_seed) 
	{
		this.rnd = new Random();
	}
	
	public double getRandomValue(double min, double max) 
	{
		double value = 0.0;
		value = this.rnd.nextDouble()*(max - min) + min;
		return value;
	}
	
	public int getRandomIndex(int size)
	{
		int index = -1;
		index = (int) Math.round(this.rnd.nextDouble() * size);
		
		if (index == size)
			index--;
		if (index == -1)
			index = 0;
		
		return index;
	}	
}
