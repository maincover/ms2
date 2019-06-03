class Df_timegrid:
    """put df and time_grid as class value then saved to pickle"""

    def __init__(self, df, relative_time_grid):
        self.df = df;
        self.time_grid_relative=relative_time_grid;
        #print df
    
    def checkLength(self):
        if self.df is not None and self.time_grid_relative is not None:
            print "length of df is: ", len(self.df);
            print "length of relative_time_grid is: ",len(self.time_grid_relative);