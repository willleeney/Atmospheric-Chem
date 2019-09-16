######################################### PLOT IN WEBGL VP ########################################################
ryb = cl.scales['3']['qual']['Set1']
legend_key = ['Original','O', 'OO', 'C=O', 'C(=O)O', '#N', 'N', 'Cl', 'Br', 'I', 'OCOO', '(=O)(OCOO)', 'OOCOO', 'O(mid)', 'O(all)']

print('What two functional groups do you want to compare?')
for i in range(len(legend_key)):
    print(f'{i} : {legend_key[i]}')
comp1 = int(input('Please choose the first group: '))
comp2 = int(input('Please choose the second group: '))


ryb = cl.scales['3']['qual']['Set1']

range_id = [0, 1175] # what range of molecules [0,1175 is all]

example = all_compounds_rate[all_compounds_rate['ID'] == comp1] #addition[(addition['Colour'] == choice[0]) & (addition['Parent']>= range_id[0]) & (addition['Parent']<= range_id[1])].copy()
example1 = all_compounds_rate[all_compounds_rate['ID'] == comp2]#addition[(addition['Colour'] == choice[1]) & (addition['Parent']>= range_id[0]) & (addition['Parent']<= range_id[1])].copy()
#example = pd.merge(example,example1, how='outer')
#N = len(example)

# add each set of plots 
called = input('What would you like to save the graphs as: ' )

trace = go.Scattergl(
    x = example.loc[:, 'Parent'],
    y = example.loc[:, 'Vapour_Pressure'],
    text = example.loc[:, 'Compound_Name'],
    mode = 'markers',
    marker = dict(
        color = ryb[0]#list(example.loc[:, 'Colour'])
    ),
    name = legend_key[comp1]

    
)
trace1 = go.Scattergl(
    x = example1.loc[:, 'Parent'],
    y = example1.loc[:, 'Vapour_Pressure'],
    text = example1.loc[:, 'Compound_Name'],
    mode = 'markers',
    marker = dict(
        color = ryb[1]#list(example1.loc[:, 'Colour'])
    ),
    name = legend_key[comp2]

    
)
data = [trace, trace1]
layout = go.Layout( # plot title
    #title=go.layout.Title(
    #   text='Plot Title',
    #    xref='paper',
    #    x=0
    #),
    xaxis=go.layout.XAxis( # xaxis label
        title=go.layout.xaxis.Title(
            text='Parent Compound ID',
            font=dict(
                family='Courier New, monospace',
                size=18,
                color='#7f7f7f'
            )
        )
    ),
    yaxis=go.layout.YAxis( #y axis label
        title=go.layout.yaxis.Title(
            text='Vapour Pressure',
            font=dict(
                family='Courier New, monospace',
                size=18,
                color='#7f7f7f'
            )
        )
    )
)
fig = go.Figure(data=data, layout=layout)

py.iplot(fig, filename=f'{called}_vp') ###################


trace2 = go.Scattergl(
    x = example.loc[:, 'Parent'],
    y = example.loc[:, 'Rate_Constant'],
    text = example.loc[:, 'Compound_Name'],
    mode = 'markers',
    marker = dict(
        color = ryb[0]#list(example.loc[:, 'Colour'])
    ),
    name = legend_key[comp1]

    
)
trace3 = go.Scattergl(
    x = example1.loc[:, 'Parent'],
    y = example1.loc[:, 'Rate_Constant'],
    text = example1.loc[:, 'Compound_Name'],
    mode = 'markers',
    marker = dict(
        color = ryb[1]#list(example1.loc[:, 'Colour'])
    ),
    name = legend_key[comp2]

    
)
data1 = [trace2, trace3]
layout1 = go.Layout( # plot title
    #title=go.layout.Title(
    #   text='Plot Title',
    #    xref='paper',
    #    x=0
    #),
    xaxis=go.layout.XAxis( # xaxis label
        title=go.layout.xaxis.Title(
            text='Parent Compound ID',
            font=dict(
                family='Courier New, monospace',
                size=18,
                color='#7f7f7f'
            )
        )
    ),
    yaxis=go.layout.YAxis( #y axis label
        title=go.layout.yaxis.Title(
            text='Rate Constant',
            font=dict(
                family='Courier New, monospace',
                size=18,
                color='#7f7f7f'
            )
        )
    )
)
fig1 = go.Figure(data=data1, layout=layout1)

py.iplot(fig1, filename=f'{called}_rc')



trace4 = go.Scattergl(
    x = example.loc[:, 'Vapour_Pressure'],
    y = example.loc[:, 'Rate_Constant'],
    text = example.loc[:, 'Compound_Name'],
    mode = 'markers',
    marker = dict(
        color = ryb[0]#list(example.loc[:, 'Colour'])
    ),
    name = legend_key[comp1]

    
)
trace5 = go.Scattergl(
    x = example1.loc[:, 'Vapour_Pressure'],
    y = example1.loc[:, 'Rate_Constant'],
    text = example1.loc[:, 'Compound_Name'],
    mode = 'markers',
    marker = dict(
        color = ryb[1]#list(example1.loc[:, 'Colour'])
    ),
    name = legend_key[comp2]

    
)
data2 = [trace4, trace5]
layout2 = go.Layout( # plot title
    #title=go.layout.Title(
    #   text='Plot Title',
    #    xref='paper',
    #    x=0
    #),
    xaxis=go.layout.XAxis( # xaxis label
        title=go.layout.xaxis.Title(
            text='Vapour Pressure',
            font=dict(
                family='Courier New, monospace',
                size=18,
                color='#7f7f7f'
            )
        )
    ),
    yaxis=go.layout.YAxis( #y axis label
        title=go.layout.yaxis.Title(
            text='Rate Constant',
            font=dict(
                family='Courier New, monospace',
                size=18,
                color='#7f7f7f'
            )
        )
    )
)
fig2 = go.Figure(data=data2, layout=layout2)

py.iplot(fig2, filename=f'{called}_vp_rc')