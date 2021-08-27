function [XTrain1,XTrain2,XTest1,XTest2] = divideTrainTest(data1,data2,indClass1,indClass2,idxFold)
%divide data into training and testing groups

testIndClass1 = (indClass1 == idxFold);
testIndClass2 = (indClass2 == idxFold);

trainIndClass1 = ~testIndClass1;
trainIndClass2 = ~testIndClass2;

if ndims(data1)>2
    
    % Get Trainingdata
    
    XTrain1 = data1(:,:,trainIndClass1);
    XTrain2 = data2(:,:,trainIndClass2);
    
    % Get Testdata
    
    XTest1 = data1(:,:,testIndClass1);
    XTest2 = data2(:,:,testIndClass2);
    
else
    
    % Get Trainingdata
    
    XTrain1 = data1(trainIndClass1);
    XTrain2 = data2(trainIndClass2);
    
    
    % Get Testdata
    
    XTest1 = data1(testIndClass1);
    XTest2 = data2(testIndClass2);
    
end



end

