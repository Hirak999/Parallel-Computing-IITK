#!/usr/bin/env python
# coding: utf-8

# In[7]:


import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import math








df=pd.read_csv('output1.txt',delimiter=' ',names=['D','P','ppn','mode','time'])


# ## Adding column (P, ppn)

# In[3]:


df["(P, ppn)"] = list(map(lambda x, y: ("(" + x + ", " + y + ")"), map(str, df["P"]), map(str, df["ppn"])))




# ## Creating -Log(time)

# In[5]:


df['2^time']=df['time']


# In[8]:


df['2^time']=df['2^time'].apply(lambda x:math.pow(2,x))


# ## Plotting against -Log(time)

# In[9]:


fig, ax = plt.subplots(1, 3, figsize=(20, 4))

for i in range(len(df['D'].unique())):
    df1=df[df['D']==df['D'].unique()[i]]
    sns.barplot(x="(P, ppn)", y="2^time", data=df1,hue='mode',ax=ax[i])
    ax[i].title.set_text('D= '+str(df['D'].unique()[i]))
#plt.show()

plt.savefig("Plot for MPI_BCAST")

















# ## Importing Files

# In[2]:


df=pd.read_csv('output2.txt',delimiter=' ',names=['D','P','ppn','mode','time'])


# ## Adding column (P, ppn)

# In[3]:


df["(P, ppn)"] = list(map(lambda x, y: ("(" + x + ", " + y + ")"), map(str, df["P"]), map(str, df["ppn"])))




# ## Creating -Log(time)

# In[5]:


df['2^time']=df['time']


# In[8]:


df['2^time']=df['2^time'].apply(lambda x:math.pow(2,x))


# ## Plotting against -Log(time)

# In[9]:


fig, ax = plt.subplots(1, 3, figsize=(20, 4))

for i in range(len(df['D'].unique())):
    df1=df[df['D']==df['D'].unique()[i]]
    sns.barplot(x="(P, ppn)", y="2^time", data=df1,hue='mode',ax=ax[i])
    ax[i].title.set_text('D= '+str(df['D'].unique()[i]))
#plt.show()
plt.savefig("Plot for MPI_REDUCE")























df=pd.read_csv('output3.txt',delimiter=' ',names=['D','P','ppn','mode','time'])


# ## Adding column (P, ppn)

# In[3]:


df["(P, ppn)"] = list(map(lambda x, y: ("(" + x + ", " + y + ")"), map(str, df["P"]), map(str, df["ppn"])))




# ## Creating -Log(time)

# In[5]:


df['2^time']=df['time']


# In[8]:


df['2^time']=df['2^time'].apply(lambda x:math.pow(2,x))


# ## Plotting against -Log(time)

# In[9]:


fig, ax = plt.subplots(1, 3, figsize=(20, 4))

for i in range(len(df['D'].unique())):
    df1=df[df['D']==df['D'].unique()[i]]
    sns.barplot(x="(P, ppn)", y="2^time", data=df1,hue='mode',ax=ax[i])
    ax[i].title.set_text('D= '+str(df['D'].unique()[i]))
#plt.show()
plt.savefig("Plot for MPI_GATHER")












df=pd.read_csv('output4.txt',delimiter=' ',names=['D','P','ppn','mode','time'])


# ## Adding column (P, ppn)

# In[3]:


df["(P, ppn)"] = list(map(lambda x, y: ("(" + x + ", " + y + ")"), map(str, df["P"]), map(str, df["ppn"])))




# ## Creating -Log(time)

# In[5]:


df['time']=df['time']


# In[8]:


df['time']=df['time'].apply(lambda x:x)



fig, ax = plt.subplots(1, 3, figsize=(20, 4))

for i in range(len(df['D'].unique())):
    df1=df[df['D']==df['D'].unique()[i]]
    sns.barplot(x="(P, ppn)", y="time", data=df1,hue='mode',ax=ax[i])
    ax[i].title.set_text('D= '+str(df['D'].unique()[i]))
#plt.show()
plt.savefig("Plot for Alltoallv");
