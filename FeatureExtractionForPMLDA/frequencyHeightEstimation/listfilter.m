function  OutputList =  listfilter(InputList, ext);
% version 0001 3-21-2005
% this programm takes a cellarray list of file names 
% and remvoes from the list those files whose extention 
% does not match the  ext string. the output is a cellarray.

% sample use  OutputList = listfilter(InputList, 'jpg') places
% only those files with .jpg extention in the OutputList

count = 0;
for index = 1:length(InputList)
      item = char(InputList(index));
      if strcmp(item(end-2:end), ext) == 1
          count = count + 1;
          indexarray(count) = index;
      end
end
for index = 1:length(indexarray)
    OutputList(index) = InputList(indexarray(index));
end

OutputList = OutputList';