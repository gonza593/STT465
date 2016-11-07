## Bonus HW (up to 10 points towards the exam)

Today I check the [Real Clear Politics](http://www.realclearpolitics.com/epolls/2016/president/2016_elections_electoral_college_map.html) 
electoral map and merge their data with that of the [electoral votes per state](https://www.archives.gov/federal-register/electoral-college/allocation.html).

The file [electoral_votes.txt](https://github.com/gdlc/STT465/blob/master/electoral_votes.txt) has the number of electoral votes per state, whether it is a state that will be almost surely carried by
Donald Trump (T) or Hillary Clinton (H) or it is currently viewed by RCP as a Toss-up (TU). For those that are TU the last
column gives the current estimate (average of several polls) of the % likely voters that manifest they will vote for Trump among 
the sum of likely voters that manifest they will vote for T or H.

Finally, although I have no real idea of the sample size used to estimate the % of likely Trump/Hillary supporters, I added a column
where I put a theoretical sample size considering that regular polls in the US include ~1000 participants and that the average reported by
RCP usually include 4-6 independent polls.

Using the data provided 

   (1) Implement a Monte Carlo algorithm for estimating the probability that Hillary Clinton will win the election. 
   Explain the underlying model and your algorithm.
   
   (2) Report the posterior distribution of the probability that Hillary Clinton will win the election
   
   (3) What is the estimated probability that Donald Trump will win the election?

Submit your HW at my e-mail gustavoc@msu.edu.

**Submit beefore Nov. 8th, 8pm  EST time**: get up to 15 points towards the miterm.

**Submit between Nov. 8th 8pm and Saturday Nov. 12th**: get up to 10 points towards the midterm.

In any case be ready to explain your approach/answeres in class.

