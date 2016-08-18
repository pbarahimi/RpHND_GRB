package RpHND_GRB.model;

import java.security.InvalidAlgorithmParameterException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Set;


public class BinaryTree {
	
	public static int getLevel(int index){
		return (int) Math.floor((Math.log(index+1)/Math.log(2)));
	}
	
	public static int getLeftChild(int index, int depth) {
		if (Math.floor(Math.log(index+1)/Math.log(2))>depth-1){
			throw new NoSuchElementException("Node " + index + " with depth " + depth + " has no children");
		}
		return 2*index + 1;
	}
	
	public static int getRightChild(int index, int depth) {
		if (Math.floor(Math.log(index+1)/Math.log(2))>depth-1){
			throw new NoSuchElementException("Node " + index + " with depth " + depth + " has no children");
		}
		return 2*index + 2;
	}
	
	public static List<Integer> getLeftChildren(int index, int depth) {
		if (getLevel(index)>depth-1){
			throw new NoSuchElementException("Node " + index + " with depth " + depth + " has no children");
		}
		int leftchild = (2*index)+1;
		
		List<Integer> result = new ArrayList<Integer>();
		List<Integer> unvisitedList = new ArrayList<Integer>();
		result.add(leftchild);
		unvisitedList.add(leftchild);
		int currentNode = unvisitedList.get(0);
		unvisitedList.remove(0);
		
		while (getLevel(currentNode)<depth){
			
			int newLeftChild = (2*currentNode)+1;
			int newRightChild = (2*currentNode)+2;
			
			result.add(newLeftChild);
			result.add(newRightChild);
			unvisitedList.add(newLeftChild);
			unvisitedList.add(newRightChild);
			currentNode = unvisitedList.get(0);
			
			unvisitedList.remove(0);
		}
		return result;				
	}
	
	public static List<Integer> getRightChildren(int index, int depth) {
		if (getLevel(index)>depth-1){
			throw new NoSuchElementException("Node " + index + " with depth " + depth + " has no children");
		}
		int rightchild = (2*index)+2;
		
		List<Integer> result = new ArrayList<Integer>();
		List<Integer> unvisitedList = new ArrayList<Integer>();
		result.add(rightchild);
		unvisitedList.add(rightchild);
		
		int currentNode = unvisitedList.get(0);
		unvisitedList.remove(0);
		
		while (getLevel(currentNode)<depth){
			int newLeftChild = (2*currentNode)+1;
			int newRightChild = (2*currentNode)+2;
			
			result.add(newLeftChild);
			result.add(newRightChild);
			unvisitedList.add(newLeftChild);
			unvisitedList.add(newRightChild);
			currentNode = unvisitedList.get(0);
			
			unvisitedList.remove(0);
		}
		return result;				
	}

	public static List<Integer> getParents(int index){
		List<Integer> result = new ArrayList<Integer>();		
		while (index>0){
			int parent = (int) (Math.floor((index-1)/2));
			result.add(parent);
			index = parent;
		}
		return result;
	}

	public static List<Integer> getLeafNodes(int root, int depth){
		List<Integer> result = new ArrayList<Integer>();
		List<Integer> temp = new ArrayList<Integer>();
		Set<Integer> toBeRemoved = new HashSet<Integer>();
		result.add(root);
		
		while(depth>0){
			for (Integer i:result){
				int leftChild = 2*i+1;
				int rightChild = 2*i+2;
				temp.add(leftChild);
				temp.add(rightChild);
				toBeRemoved.add(i);
			}
			result.removeAll(toBeRemoved);
			result.addAll(temp);
			temp.clear();
			toBeRemoved.clear();			
			depth--;			
		}
		return result;
	}

	public static boolean isLeftChild(int parent, int child) throws InvalidAlgorithmParameterException{
		int depth = getLevel(child) - getLevel(parent);
		List<Integer> parentsLeafNodes = getLeafNodes(parent, depth);
		if (!parentsLeafNodes.contains(child)){
			throw new InvalidAlgorithmParameterException(child + " is not a leaf node for " + parent);
		}
		
		if (child < (double) (parentsLeafNodes.get(0) + parentsLeafNodes.get(parentsLeafNodes.size()-1))/2){
			return true;
		}
		else{
			return false;
		}
	}
}
