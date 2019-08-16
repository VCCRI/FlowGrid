import numpy as np
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import MinMaxScaler
from copy import deepcopy


class div_bin():
    def __init__(self,X,bin_n):
        self.data=X
        self.bin_n=bin_n

    def DividingOneD(self,X,d2,index=None):
        if index is None:
            index=np.arange(X.shape[0])
        x=X[index,d2]
        ############################################
        #print(self.bin_n)
        bin_list=[None for i in range(self.bin_n)]
        ############################################
        #print(bin_list)
        temp_index=np.arange(x.shape[0])
        st,end=np.min(x),np.max(x)
        for i in range(st,end+1):
            cur_index=np.equal(x[temp_index],i)
            tem_index=temp_index[cur_index]
            ############################################
            #print(tem_index)
            #print(index)
            if tem_index.shape[0] >0:
                bin_list[i]=index[tem_index]
                temp_index=temp_index[np.logical_not(cur_index)]
        ############################################        
        #print(bin_list)
        return bin_list

    def combine(self,index_list,d2):
        X=self.data
        new_index_dict={}
        for index,value in enumerate(index_list):
            if value is not None:
                th_index_list=self.DividingOneD(X,d2,value)
                for k,i in enumerate(th_index_list):
                    if i is not None:
                        new_index_dict[(index,k)]=i
        return new_index_dict

    def dividing_bins(self):
        X=self.data
        d=X.shape[1]
        th_index_list=self.DividingOneD(X,0)
        new_index_dict=self.combine(th_index_list,1)
        if d>2:
            next_d={}
            for key,value in new_index_dict.items():
                temp_dict={key:value}
                for cur_d in range(2,d):
                    new_d={}
                    for cur_key,cur_index in temp_dict.items():
                        th_index_list=self.DividingOneD(X,cur_d,cur_index)
                        for k,i in enumerate(th_index_list):
                            if i is not None:
                                new_d[cur_key+(k,)]=i
                    if cur_d!=d-1:
                        temp_dict=deepcopy(new_d)
                for k,i in new_d.items():
                    next_d[k]=i
            return next_d
        else:
            return new_index_dict

    def unique(self):
        idd=self.dividing_bins()
        id_list=list(idd.keys())
        counts=np.zeros(len(id_list),dtype=int)
        unique_index=np.zeros(self.data.shape[0],dtype=int)
        for k,i in enumerate(id_list):
            counts[k]=idd[i].shape[0]
            unique_index[idd[i]]=k
        return np.array(id_list),unique_index,counts

class FlowGrid():
    def __init__(self,original_data,MinDenB=3,bin_n=14,eps=1.5, MinDenC=40):
        """
		Here is the initial step of building the object of FlowGrid.
		MinDenB, max_value, bin_n, eps and MinDenC are set in this step.
		1. MinDenB is the minimun density for high density bin, which 
		could be changed by user-input.
		2. max_value stands for the range of normalized data, ranging 
		from zero to five.
		3. bin_n stands for the number of equal sized bins for each dimension.
		4. eps is the maximun distance of connected bins.
		5. MinDenC stands for the minimun collective density for core bins.
		After setting the five parameter, _initial_ is used to normalize
		the original data to normalized data ranging from zero to five and 
		calculate the bandwith of each bin based on the setting of bin_n (
		number of bin).
		"""
        self.original_data=np.array(original_data)
        self.MinDenB=MinDenB
        self.bin_n=bin_n
        self.eps=eps
        self.MinDenC=MinDenC
        self.n,self.d=self.original_data.shape

    def dividing_bins(self):
        """
        The function is to divid the data into equal sized bins.
        It starts with normalizing data ranging from 0 to bin_n-10**(-6).
        The integer part of each value could represent the coordinate of data 
        so we first change the type of data into int and np.unique could 
        get the distint set of coordinate with sample id and number of cell 
        locating in the bin.
        """
        scaler = MinMaxScaler(feature_range=(0, self.bin_n-10**(-7)), copy=False)
        scaler.fit(self.original_data)
        X=scaler.transform(self.original_data).astype(u"int8")
        unique_array,unique_index,counts=div_bin(X,self.bin_n).unique()

        if np.log(self.n)/np.log(10)>4 and self.d<11 and self.bin_n<11:
            unique_array,unique_index,counts=div_bin(X,self.bin_n).unique()
        else:
            unique_array,unique_index,counts=np.unique(X,return_inverse=True,return_counts=True, axis=0)
        return unique_array,unique_index, counts

    def density_query(self,unique_array,counts,nn_model):
        """
        The density query function is to determine whether a bin is a core
        or not. 
        It starts with filtering out the bins whose density is lower than
        MinDenB, followed by near neighbor search using radius_neighbors from
        Sklearn.
        The next step is to calculate the collective density if the number of 
        samples located in the bin is larger than 85% of its connected bins.
        """
        self.bins_number=unique_array.shape[0]
        tf_array=np.greater(counts,self.MinDenB)
        index_array=np.arange(self.bins_number)
        check_index=index_array[tf_array]
        check_nn=unique_array[tf_array]
        filterd_size=check_nn.shape[0]
        neighborhoods = nn_model.radius_neighbors(check_nn,self.eps,return_distance=False)
        n_neighbors =np.zeros((filterd_size))
        for k,neighbors in enumerate(neighborhoods):
            nn_list=counts[neighbors]
            key_n=counts[check_index[k]]
            if key_n>=np.percentile(nn_list,85):
                n_neighbors[k]=np.sum(nn_list)
        core_non=np.where(n_neighbors>=self.MinDenC)[0]
        core_or_non=np.zeros(self.bins_number,dtype=bool)
        core_bin_index=check_index[core_non]
        core_or_non[core_bin_index]=True
        query_d={}
        for core in core_non:
            query_d[check_index[core]]=neighborhoods[core] 
        return query_d,core_or_non

    def bfs(self,query_d,core_non):
        """
        bfs (breadth first search) is to group the core bins with
        their connected bins.
        The initial setting is assigning all bin label as -1 which
        stands for noise and building a set object for a queue.
        In the first while loop, if the core bin has not been labelled,
        it will be pushed into the queue, while if it is labelled,
        filter function is applied to remove all the labelled bins
        in bin_list.
        The second level while loop is to connect the core bins with 
        their connected bins. cur_bin is poped from queue. If it is 
        not be labelled, it will be label by index and if it is core 
        bin, all non-labelled bin is put into queue for the next iteration.
        """
        bin_labels= np.zeros(self.bins_number)-1
        filter_f=lambda x: bin_labels[x]==-1
        index=0
        queue=set()
        bin_list=list(query_d.keys())
        while bin_list:
            core_bin=bin_list.pop()
            if bin_labels[core_bin]==-1:
                queue.add(core_bin)
                index+=1
                while len(queue)>0:
                    cur_bin=queue.pop()
                    if bin_labels[cur_bin]==-1:
                        bin_labels[cur_bin]=index
                        if core_non[cur_bin]:
                            queue.update(list(filter(filter_f,query_d[cur_bin])))
            else:
                bin_list=list(filter(filter_f,bin_list))
        return bin_labels

    def density_scan(self,unique_array, counts):
        """
        The function is to group core bins and their connected bins.
        It starts with querying neighbor bins by applying NearestNeighbors
        function in sklearn. The function applies kd-tree to reduce the time
        complixity from O(n^2) to O(nlogn).
        And them, the density_query is to determine the core bins and non-core
        bins while bfs(breath first search) function is to grouping the connected 
        bins.
        """
        nn_model = NearestNeighbors()
        nn_model.fit(unique_array)
        neighbors_d,core_non=self.density_query(unique_array,counts,nn_model)
        bin_labels=self.bfs(neighbors_d,core_non)
        return bin_labels

    def clustering(self):
        #clustering function is the main function calling the above functions.
        #It starts with dividing data into equal sized bins.
        unique_array,unique_index, counts=self.dividing_bins()
        #It is followed by searching core bins and grouping them with their near bins.
        bin_labels=self.density_scan(unique_array, counts)
        #Labeling sample by the label of corresponding bins.
        label_array=bin_labels[unique_index]
        return label_array
