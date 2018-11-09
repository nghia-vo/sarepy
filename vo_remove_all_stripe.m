% #============================================================================
% # Description: Matlab implementation of stripe artifact removal methods, 
% # Nghia T. Vo, Robert C. Atwood, and Michael Drakopoulos, "Superior
% # techniques for eliminating ring artifacts in X-ray micro-tomography," Optics
% # Express 26, 28396-28412 (2018).
% # https://www.osapublishing.org/oe/fulltext.cfm?uri=oe-26-22-28396&id=399265
% #============================================================================
% Translation to Matlab R2018a: 
% Translator: Daniel S. Hussey
% email: daniel.hussey@nist.gov
% Date: 08 NOV 2018
% Translation only of algorithms 3,4,5,6
% sinograms have angular coordinates along dimension 1

function [sinogram] = vo_remove_all_stripe(sinogram, snr, la_size, sm_size, filterstr)
%     Remove all types of stripe artifacts by combining algorithm 6, 5, and 3.
%     Angular direction is along the axis 0.
%     ---------
%     Parameters: - sinogram: 2D array.
%                 - snr: ratio used to discriminate between useful
%                     information and noise
%                 - la_size: size of the median filter to remove
%                     large stripes.
%                 - sm_size: size of the median filter to remove
%                     small-to-medium stripes.
%                 - filters - a string, if set to 'all' run algorithms
%                     3, 6, and  5 otherwise, just run algorithms 3 and 5
%     ---------
%     Return:     - stripe-removed sinogram.
    % Run algorithms 3, 5 and/or 6, algortihm 4 is called by 5 and 6
    sinogram = remove_stripe_based_sorting(sinogram, sm_size);
    if contains(filterstr,'all','IgnoreCase',true)
        sinogram = remove_unresponsive_or_fluctuating_stripe(sinogram, snr, la_size);
    else
        sinogram = remove_large_stripe(sinogram, snr, la_size);
    end
    
end

function [sino_corrected, sino_sort, indices, filt_sino_sort] = remove_stripe_based_sorting(sinogram, sm_size)
%     Algorithm 3 in the paper. Remove stripes using the sorting technique.
%     Work particularly well for removing partial stripes.
%     Angular direction is along the axis 0
%     ---------
%     Parameters: - sinogram: 2D array. 
%                 - size: size of the median filter.
%     ---------
%     Return:     - stripe-removed sinogram.
    [nrow, ncol] = size(sinogram);
    [sino_sort, indices] = sort(sinogram, 1);
    filt_sino_sort = medfilt2(sino_sort,[sm_size,sm_size], 'symmetric');
    for jx = 1:ncol
        sino_corrected(indices(:,jx),jx) = filt_sino_sort(1:nrow,jx);
    end
end
 
    function [listmask] = detect_stripe(listdata, snr)
%     Algorithm 4 in the paper. Used to locate stripe positions.
%     ---------
%     Parameters: - listdata: 1D normalized array.
%                 - snr: ratio used to discriminate between useful
%                     information and noise
%     ---------
%     Return:     - 1D binary mask.
    numdata = numel(listdata);
    listsorted = sort(listdata);
    xlist = (linspace(1,numdata,numdata));
    ndrop = floor(0.25*numdata);
    fitresults = fit(double(transpose(xlist(ndrop:numdata-ndrop-1))),double(transpose(listsorted(ndrop:numdata-ndrop-1))),'poly1');
    vals = single(coeffvalues(fitresults));
    slope = vals(1); intercept = vals(2);
    numt1 = intercept + slope * xlist(numdata);
    noiselevel = abs(numt1 - intercept);
    val1 = abs(listsorted(1) - intercept) / noiselevel;
    val2 = abs(listsorted(numdata) - numt1) / noiselevel;
    listmask = zeros(1,numdata);
    if (val1 >= snr)
        upper_thresh = intercept + noiselevel * snr * 0.5;
        listmask(listdata > upper_thresh) = 1.0;
    end
    if (val2 >= snr)
        lower_thresh = numt1 - noiselevel * snr * 0.5;
        listmask(listdata <= lower_thresh) = 1.0;
    end
end

function [sinogram] = remove_large_stripe(sinogram, snr, la_size)
%     Algorithm 5 in the paper. Remove large stripes by: locating stripes,
%     normalizing to remove full stripes, using the sorting technique to
%     remove partial stripes.
%     Angular direction is along the axis 0.
%     ---------
%     Parameters: - sinogram: 2D array.
%                 - snr: ratio used to discriminate between useful
%                     information and noise
%                 - size: size of the median filter.
%     ---------
%     Return:     - stripe-removed sinogram.
    badpixelratio = 0.05;
    [nrow, ncol] = size(sinogram);
    [sino_corrected, sinosorted, indices, sinosmoothed] = remove_stripe_based_sorting(sinogram, la_size);
    ndrop = floor(badpixelratio * nrow);
    list1 = mean(sinosorted(ndrop:nrow - ndrop,:), 1);
    list2 = mean(sinosmoothed(ndrop:nrow - ndrop,:),1);
    listfact = list1 ./ list2;
    listmask = detect_stripe(listfact, snr);
    listmask = imdilate(listmask,[1,1,1]);
    matfact = repmat(listfact,[nrow,1]);
    sinogram = sinogram ./matfact;
    [sinogram_step4, ~] = remove_stripe_based_sorting(sinogram, la_size);
    sinogram(:, listmask > 0.0) = sinogram_step4(:, listmask > 0.0);
end

function [sinogram] = remove_unresponsive_or_fluctuating_stripe(sinogram, snr, la_size)
%     Algorithm 6 in the paper. Remove unresponsive or fluctuating stripes by:
%     locating stripes, correcting by interpolation.
%     Angular direction is along the axis 0.
%     ---------
%     Parameters: - sinogram: 2D array.
%                 - snr: ratio used to discriminate between useful
%                     information and noise
%                 - size: size of the median filter.
%     ---------
%     Return:     - stripe-removed sinogram.
    [nrow, ncol] = size(sinogram);
    filter1d = ones(1,11)./11;
    sinosmoothed = filter2(filter1d,padarray(sinogram,[0,5],'replicate','both' ),'valid');
    listdiff = sum(abs(sinogram - sinosmoothed),1);
    nmean = mean(listdiff);
    listdiffbck = medfilt1(listdiff, la_size);
    listdiffbck(listdiffbck == 0.0) = nmean;
    listfact = listdiff ./ listdiffbck;
    listmask = detect_stripe(listfact, snr);
    listmask = imdilate(listmask,[1,1,1]);
    listmask(1:3) = 0; listmask(nrow-2:nrow) = 0;
    xind = transpose(linspace(1,ncol,ncol));
    yind = transpose(linspace(1,nrow,nrow));
    listx = xind(listmask < 1.0);
    matz = sinogram(:,listx);
    [X, Y] = meshgrid(listx, yind);
    [Xq, Yq] = meshgrid(xind, yind);
    finter = interp2(X, Y, matz, Xq, Yq, 'linear');
    if numel(listmask(listmask > 0)) > 0
        sinogram(:, xind(listmask > 0)) = finter(:, xind(listmask > 0));
    end
    %      Use algorithm 5 to remove residual stripes
    sinogram = remove_large_stripe(sinogram, snr, la_size); 

 end

