module Gating

    using Dash, DataFrames, JSON, PlotlyJS, FlowCytometry

    function manualGating(fcs::FlowCytometryExperiment)

        gates = Dict{String,Any}()

        app = dash()

        dropdown_options = [
            Dict("label" => i, "value" => i) for i in fcs.channels
        ]

        function generate_table(gates)
            return [
                html_thead(html_tr([html_th("Name"),html_th("ChannelX"),html_th("ChannelY"),html_th("Type")])),
                html_tbody([html_tr([html_td(r), html_td(gates[r]["channels"][1]), html_td(gates[r]["channels"][2]), html_td(gates[r]["type"])]) for r = keys(gates)]),
            ]
        end

        p1 = Plot(fcs[fcs.channels[1]], fcs[fcs.channels[2]], marker_size=8)

        app.layout = html_div() do
            html_h1("Gating",
                style = Dict("color" => "#00000", "textAlign" => "center")),
            html_div("Dash: A web application framework for your data."),
            html_div(style = Dict("columnCount" => 2)) do
                dcc_graph(
                    id = "scatter",
                    figure = p1
                ),
                html_div(html_th("Visualizing gates")),
                html_div(style = Dict("columnCount" => 2)) do
                    html_label("Gate X"),
                    dcc_dropdown(
                        id = "inputX",
                        options = dropdown_options,
                        value = names(d)[1],
                        multi = false,
                    ),
                    html_label("Gate Y"),
                    dcc_dropdown(
                        id = "inputY",
                        options = dropdown_options,
                        value = names(d)[1],
                        multi = false,
                    )
                end,
                html_div(html_th("Add gate")),
                html_div(style = Dict("columnCount" => 2)) do
                    dcc_input(id = "gate-input-name", value = "Add name of gate here...", type = "text"),
                    html_button(id = "submit-gate", children = "add gate", n_clicks = 0),
                    html_div(id="added-gate-check","")
                end,
                html_label(html_th("Remove gate")),
                html_div(style = Dict("columnCount" => 2)) do
                    dcc_dropdown(
                        id = "gate-remove-name",
                        options = keys(gates),
                        value = "",
                        multi = false,
                    ),
                    html_button(id = "delete-gate", children = "remove gate", n_clicks = 0),
                    html_div(id="remove-gate-check","")
                end,
                html_div(id="relayout-data"),
                html_table(
                    id = "gates",
                    [html_thead(html_tr(["Name","ChannelX","ChannelY"]))]
                )
            end
        end

        callback!(app, Output("scatter", "figure"), Input("inputX", "value"), Input("inputY", "value")) do x,y
            fig = plot(
                scatter(x=fcs[x], y=fcs[y], mode="markers", marker_size=8, yaxis_title="Hola"),
                Layout(xaxis_title=x,yaxis_title=y,height=800,width=800)
            )
            fig
        end

        callback!(app, Output("added-gate-check","children"), Input("submit-gate", "n_clicks"), State("gate-input-name","value"), State("inputX", "value"), State("inputY", "value"), State("scatter", "selectedData")) do button,name,x,y,box
            if name in keys(gates)
                return "Gate name already exisiting in gates."
            elseif "Add name of gate here..." == name
                return ""
            else
                dict = Dict(box)
                gates[name] = Dict{String,Any}()
                gates[name]["channels"] = (x,y)
                gates[name]["points"] = [i[:pointNumber] for i in box[:points]]
                if "range" in keys(box)
                    gates[name]["type"] = "box"
                    gates[name]["boundary"] = Dict{String,Vector}()
                    gates[name]["boundary"]["x"] = box[:range][:x]
                    gates[name]["boundary"]["y"] = box[:range][:y]
                else
                    gates[name]["type"] = "lasso"
                    gates[name]["boundary"] = Dict{String,Vector}()
                    gates[name]["boundary"]["x"] = box[:lassoPoints][:x]
                    gates[name]["boundary"]["y"] = box[:lassoPoints][:y]
                end
                return "Gate added!"
            end
        end

        callback!(app, Output("remove-gate-check","children"), Input("delete-gate", "n_clicks"), State("gate-input-name","value")) do button,name
            delete!(gates,name)
            return ""
        end

        callback!(app, Output("gates","children"), Output("gate-remove-name","options"), Input("submit-gate", "n_clicks"), Input("delete-gate", "n_clicks")) do addbutton,delbutton
            return generate_table(gates), [Dict("label" => i, "value" => i) for i in keys(gates)]
        end


        run_server(app, "0.0.0.0", debug=true)

        return

    end

end