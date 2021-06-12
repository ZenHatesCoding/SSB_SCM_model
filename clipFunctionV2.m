function [clippedSignal]= clipFunctionV2(unclippedSignal, clipProb, numLevels)

    clipMinRatio = clipProb;
    clipMaxRatio = 1 - clipMinRatio;
    sigDist = hist(unclippedSignal,numLevels);
    maxVal = max(unclippedSignal);
    minVal = min(unclippedSignal);
    % maxValHit= max(sigDist);
    % lowHitRangesToClip = find(sigDist < (maxValHit/500));
    cdfSig = cumsum(sigDist);
    cdfSig = cdfSig/cdfSig(end);
    sigRanges = minVal: (maxVal - minVal)/((numLevels-1)): maxVal;
    clipNegIndex = cdfSig < clipMinRatio;
    clipPosindex = cdfSig > clipMaxRatio;
    clipFlagON = 1;
    if(sum(clipNegIndex) == 0) || sum(clipPosindex) == 0
        clipFlagON = 0;
    else
    clipValueNeg = max(sigRanges(clipNegIndex));
    clipValuePos = min(sigRanges(clipPosindex));
    clipValue = max(abs(clipValueNeg), clipValuePos);
    end


    clippedSignal = unclippedSignal;
    if clipFlagON == 1
        clippedSignal(clippedSignal>clipValue) = clipValue;
        clippedSignal(clippedSignal<-clipValue) = -clipValue;

    end
end