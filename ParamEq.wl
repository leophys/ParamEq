(* ::Package:: *)

(* ::Title:: *)
(*A parametric equalizer*)


(* ::Chapter:: *)
(*Author: Leonardo Barcaroli*)


(* ::Section:: *)
(*Constants*)


buffLen=1024;


SampleFreq=44100;


overlap=32;


overlap2=128;


smoothingbuffer=128;


(*freqCode=AssociationThread[Table[i,{i,buffLen}],Table[N[SampleFreq*Log[2,(i+buffLen)/buffLen]],{i,buffLen}]];*)


freqCode=AssociationThread[Table[i,{i,buffLen}],Table[N[(SampleFreq*i)/buffLen],{i,buffLen}]];


window=Table[Window[i,buffLen,smoothingbuffer],{i,buffLen}];


(* ::Section:: *)
(*Functions*)


Window[i_,buffer_,overlap_,pad_:0]:=
1/(1+Exp[-16/(overlap-pad) (i-(overlap-pad)/2)]) 1/(1+Exp[16/(overlap-pad) (i-(buffer-(overlap-pad)/2))]);


Filter[\[Omega]_,Q_,\[Omega]0_,A_]:=
N[1.+Sinh[A]((1+2 Q) ( \[Omega] \[Omega]0)/(Q (\[Omega]^2+(\[Omega] \[Omega]0)/Q+\[Omega]0^2)))^Q]


HighShelf[\[Omega]_,Q_,\[Omega]0_,floor_]:=
If[Q==0||floor==0,1,N[1/2 (1-(1-10^floor)Tanh[Q(\[Omega]-\[Omega]0)])]]


LowShelf[\[Omega]_,Q_,\[Omega]0_,floor_]:=
If[Q==0||floor==0,1,N[1/2 (1+(1-10^floor)Tanh[Q(\[Omega]-\[Omega]0)])]]


NormalizeBuffered[waveform_,bufferPow_,overlap_: 0]:=Partition[waveform,2^bufferPow,2^bufferPow-overlap,1,0]


SmoothingWindow[buffPow_,overlap_:32,pad_:0]:=
Block[{buffer=2^buffPow},
Table[
N[1/(1+Exp[-16/(overlap-pad) (i-(overlap-pad)/2)]) 1/(1+Exp[16/(overlap-pad) (i-(buffer-(overlap-pad)/2))])]
,{i,buffer}]
]


FFT::buffoverflow="The buffer `1` is longer than the sample.";
FFT[waveform_,bufferPow_,smoothingbuffer_,overlap_:0]:=
Catch[
Block[{L,buffer=2^bufferPow,tmpwaveform,cycles,fft,padding,window},
L=Length[waveform];
If[2^bufferPow>L,
Message[FFT::buffoverflow,buffer];
Abort[];
];
tmpwaveform=NormalizeBuffered[waveform,bufferPow,overlap];
cycles=Length[tmpwaveform];
window=SmoothingWindow[bufferPow,smoothingbuffer];
fft=Map[Fourier[# window]&,Take[tmpwaveform,cycles]];
Return[{fft,L,overlap}];
]
]


InverseFFT[fourier_,bufferPow_]:=
Catch[
Block[{kfunc,waveform,buffer=2^bufferPow,Lwave,overlap,cycles},
{kfunc,Lwave,overlap}=fourier;
cycles=Length[kfunc];
waveform=Map[Re[InverseFourier[#]]&,Take[kfunc,cycles]];
(*waveform=SmoothBufferEdges[waveform,bufferPow,overlap];*)
Return[(Flatten@Join@Table[waveform[[i,1;;buffer-overlap]],{i,cycles}])[[1;;Lwave]]
]
]
]


KernelFilter[filter_,params__,freqcode_,bufflen_]:=
Block[{kernel},
kernel=Table[N[filter[freqcode[i],##]]&@@params,{i,bufflen}];
Return[kernel];
]


InverseFFTK::kerlen="Fourier kernel is not of the same length as the buffer `1`";
InverseFFTK[fourier_,bufferPow_,kernel_]:=
Catch[
Block[{kfunc,waveform,buffer=2^bufferPow,Lwave,overlap,cycles},
If[
Length[kernel]!=buffer,
Message[InverseFFTK::kerlen,buffer];
Abort[];
];
{kfunc,Lwave,overlap}=fourier;
cycles=Length[kfunc];
waveform=Map[Re[InverseFourier[# kernel]]&,Take[kfunc,cycles]];
(*waveform=SmoothBufferEdges[waveform,bufferPow,overlap];*)
Return[(Flatten@Join@Table[waveform[[i,1;;buffer-overlap]],{i,cycles}])[[1;;Lwave]]
]
]
]


TrackColorN[n_Integer]:=
If[n>1,
Table[
N[{
Sqrt[2 \[Pi]]/3 PDF[NormalDistribution[2N[(i-1)/(n-1)]-1,1/3]][-1],
Sqrt[2 \[Pi]]/3 PDF[NormalDistribution[2N[(i-1)/(n-1)]-1,1/3]][0],
Sqrt[2 \[Pi]]/3 PDF[NormalDistribution[2N[(i-1)/(n-1)]-1,1/3]][1]
}]
,{i,n}],
{{1,0,0}}
]


TrackColor[n_Integer]:=
RGBColor[#]&/@TrackColorN[n]


BufferPlayer::samplelen="The waveforms are not of the same length!";

BufferPlayer[waveforms__,bufferPow_,overlap_:0,pos_]:=
Block[{len=Length[waveforms[[1]]],Nw=Length[waveforms],bufflen=2^bufferPow,tmpwaves={}},
tmpwaves=NormalizeBuffered[#,bufferPow,overlap]&/@waveforms;
ListLinePlot[
Table[
tmpwaves[[wavenum,pos]]
,{wavenum,Nw}]
,Joined->True
,PlotRange->{Full,{-1,1}}
,PlotStyle->TrackColor[Nw]
,Filling->Axis
]/;And@@Flatten@Map[Equal[Length[#],len]&,waveforms]
]


Player[waveforms__,buffer_,pos_]:=
Block[{len=Length/@waveforms,maxlen,tmpwaveforms,Nw=Length[waveforms]},
maxlen=Max[len];
If[And@@Equal[#,maxlen]&/@len,
tmpwaveforms=waveforms,
tmpwaveforms=ArrayPad[#,{0,maxlen-Length[#]},0]&/@waveforms
];
Print[Length[tmpwaveforms]];Abort[];
ListPlot[
(Part[#,pos;;pos+buffer]&/@tmpwaveforms)
,Joined->True
,PlotRange->{Full,{-1,1}}
,PlotStyle->TrackColor[Nw]
,Filling->Axis]
]



FXKernel::arglen="Argument `1` is not a list of length `2`\.1d";

FXKernel[Q__,\[Omega]0__,A__,n_Integer:1,{QH_,\[Omega]0H_,fH_},{QL_,\[Omega]0L_,fL_},overlap_,bufferPow_,freqcode_]:=
Block[{buffer=2^bufferPow,filterParams,Qfilters,highShelf,lowShelf},
If[
Or@@(Unequal[Length[#],n]&/@{Q,\[Omega]0,A}),
Message[FXKernel::arglen,MapThread[Which,Flatten[({Length[#]!=n,#})&/@{Q,\[Omega]0,A},1],0],
n
];
Abort[]
];
filterParams={Q,\[Omega]0,A};
Qfilters=
MapThread[KernelFilter[Filter,{#1,#2,#3},freqcode,buffer]&,filterParams,1];
highShelf=KernelFilter[HighShelf,{QH,\[Omega]0H,fH},freqcode,buffer];
lowShelf=KernelFilter[LowShelf,{QL,\[Omega]0L,fL},freqcode,buffer];
Return[
{Times@@Qfilters*highShelf*lowShelf,
{lowShelf,##,highShelf}&@@Qfilters
}
];
]



VisualEqualizer[fft_,len_,kernel_,filters__,freqcode_,bufferPow_,overlap_:32,smoothingbuffer_:128,pos_]:=
Block[{kFuncs,Lfilters=Length[filters]},
kFuncs=
Flat@Join[{
fft[[pos]],
fft[[pos]]kernel,
kernel
}
,filters
];
kFuncs=(MapThread[{#1,#2}&,{freqcode,#},1]&/@kFuncs);
ListLogLogPlot[
kFuncs
,Joined->True
,PlotStyle->Flatten@{Red,Blue,Green,TrackColor[Length[filters]]}
,PlotRange->{Full,{10^-5,10}}
,Filling->Table[i+3->{3},{i,Lfilters}]
]
]



(* ::Section:: *)
(*The test interface*)


ParamEqualizer[sourcewave_,buffLen_:buffLen,overlap_:overlap,smoothingbuffer_:smoothingbuffer]:=
DynamicModule[
{
Qfilters={Table[1,{i,buffLen}]},
filters,
Q1d=0,
Q2d=0,
Q3d=0,
\[Omega]1d=2,
\[Omega]2d=2.7,
\[Omega]3d=3.5,
A1d=0,
A2d=0,
A3d=0,
QHd=0,
\[Omega]Hd=3.3,
fHd=0,
QLd=0,
\[Omega]Ld=1.3,
fLd=0,
fft,
vfft,
vffttrans,
pos,
transwave=Sound@SampledSoundList[sourcewave,SampleFreq],
updatesound,
ker
},
(*transwave:=Dry[L];*)
fft=FFT[sourcewave,10,smoothingbuffer,overlap];
updatesound[ker_]:=Sound@SampledSoundList[InverseFFTK[fft,10,ker],SampleFreq];
(*transwave=updatesound;*)
Panel[
Column[{
Row[{"Buffer position",
Slider[Dynamic@pos,{1,Floor[fft[[2]]/(buffLen-overlap)]}],
Dynamic@pos,Spacer[50],Dynamic@Button["Get sound!",transwave=updatesound[Qfilters[[1]]],Method->"Queued"]}],
Panel@Column[{
Dynamic[

Qfilters=FXKernel[{10^Q1d,10^Q2d,10^Q3d},{2*10^\[Omega]1d,2*10^\[Omega]2d,2*10^\[Omega]3d},{A1d,A2d,A3d},3,{QLd,2*10^\[Omega]Ld,fLd},{QHd,2*10^\[Omega]Hd,fHd},overlap,10,freqCode];

filters=Transpose[{Values[freqCode],#}]&/@({Qfilters[[1]],##}&@@Qfilters[[2]]);
vfft=Transpose[{Values[freqCode],Abs[fft[[1,pos]]]}];
vffttrans=Transpose[{Values[freqCode],Abs[fft[[1,pos]]]*Qfilters[[1]]}];

Dynamic@ListLogLogPlot[
{filters[[1]],vfft,vffttrans}
,Joined->True
,PlotStyle->TrackColor[3]
,PlotRange->{{20,20000},{10^-10,10}}
,ImageSize->700]

]
},Center]
,
Deploy@Panel[
Column[{
Row[{
Column[{
Row[{"Q1",Slider[
Dynamic@Q1d,{0,3}]
,Dynamic@Q1d}],
Row[{"\[Omega]1",Slider[
Dynamic@\[Omega]1d,{1,4}]
,Dynamic@(2*20^\[Omega]1d)}],
Row[{"A1",Slider[
Dynamic@A1d,{-1,2}],
Dynamic@N[Sinh[A1d]]}]
}],
Spacer[50],
Column[{
Row[{"Q2",Slider[
Dynamic@Q2d,{0,3}]
,Dynamic@Q2d}],
Row[{"\[Omega]2",Slider[
Dynamic@\[Omega]2d,{1,4}]
,Dynamic@(2*20^\[Omega]2d)}],
Row[{"A2",Slider[
Dynamic@A2d,{-1,2}],
Dynamic@N[Sinh[A2d]]}]
}],
Spacer[50],
Column[{
Row[{"Q3",Slider[
Dynamic@Q3d,{0,3}]
,Dynamic@Q3d}],
Row[{"\[Omega]3",Slider[
Dynamic@\[Omega]3d,{1,4}]
,Dynamic@(2*20^\[Omega]3d)}],
Row[{"A3",Slider[
Dynamic@A3d,{-1,2}],
Dynamic@N[Sinh[A3d]]}]
}]
}],
Row[{
Column[{
Row[{"QL",Slider[
Dynamic[QLd],{-1,0}]
,Dynamic@QLd}],
Row[{"\[Omega]L",Slider[
Dynamic[\[Omega]Ld],{1,4}]
,Dynamic@(2*10^\[Omega]Ld)}],
Row[{"floor_L",Slider[
Dynamic[fLd],{0,-30}]
,Dynamic@fLd}]
}],
Spacer[100],
Column[{
Row[{"QH",Slider[
Dynamic[QHd],{-.02,0}]
,Dynamic@QHd}],
Row[{"\[Omega]H",Slider[
Dynamic[\[Omega]Hd],{1,4}]
,Dynamic@(2*10^\[Omega]Hd)}],
Row[{"floor_H",Slider[
Dynamic[fHd],{0,-30}]
,Dynamic@fHd}]
}]
}]
},Center]
],
Row[{
Dynamic@Refresh[transwave,TrackedSymbols:>{transwave}]
}]
},Center]
]
]



(* ::Section:: *)
(*A simple buffer player*)


Player[wave1_,wave2_]:=
Manipulate[
BufferPlayer[{wave1,wave2},10,overlap,i] 
,{i,1,Length[NormalizeBuffered[wave1,10,overlap] ],1}]
