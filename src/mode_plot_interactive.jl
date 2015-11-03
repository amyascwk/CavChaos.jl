# #############################################################################
# #############################################################################
#Interactive PSS plot
function plotpss_interactive(resultsdir::AbstractString;kwargs...)
    
    # #########################################################################
    #Initiate plot
    
    #Prepare plot
    fig = plt[:figure](figsize=(9,9))
    axes = fig[:gca](aspect="equal")
    
    #Plot the non-interactive version
    #and save PSS plot state
    pssarray::Array{Float64,2},residarray::Array{Int64,2} = plotpss(resultsdir;kwargs...)

    
    # #########################################################################
    #Interactions tasks
    
    #Highlight bouncepoints on PSS
    #This makes the PSS orbit from a single set of initial conditions stand out
    #from other points, and any periodicity in the orbit to be more obvious
    function pss_highlight(xdata::Float64,ydata::Float64)
        #Get appropriate entry in PSS data arrays based on position on PSS
        range::Array{Float64,1} = [plt[:axis]()...]
        thetaindex::Int64 = getbinindex(xdata,range[1],range[2],size(pssarray,2))
        chiindex::Int64 = getbinindex(ydata,range[4],range[3],size(pssarray,1))
        
        #Check that entry index is within bounds
        if thetaindex != 0 && chiindex != 0
            #Read info from data arrays and recorded results
            score::Float64 = pssarray[chiindex,thetaindex]
            resultid::Int64 = residarray[chiindex,thetaindex]
            if resultid > 0
                bouncepoints::Array{Float64,2} = readresults(resultsdir,resultid,"bouncepoints")
                
                #Plot bouncepoints on PSS as red squares with thick black borders for visibility
                plt[:plot](bouncepoints[:,1],bouncepoints[:,2],
                           ".",mec="k",mfc="w",mew=2,ms=10,axes=axes)
                plt[:axis](range)
                plt[:draw]()
            end
        end
    end
    
    #Dehighlight selected bouncepoints
    function pss_dehighlight()
        #Get lines in axes and remove them all
        lines = axes[:get_lines]()
        for line in lines
            line[:remove]()
        end
        plt[:draw]()
    end
    
    #Clear PSS plot and replot
    function pss_replot(newrange::Array{Float64,1})
        axes[:cla]()
        pssarray,residarray = plotpss(resultsdir::AbstractString;axes=axes,range=newrange)
        plt[:axis](newrange)
        plt[:draw]()
    end
    
    #Pan
    function pss_pan(dxdata::Float64,dydata::Float64)
        #Get transformed range
        range::Array{Float64,1} = [plt[:axis]()...]
        newrange::Array{Float64,1} = [range[1]-dxdata, range[2]-dxdata,
                                      range[3]-dydata, range[4]-dydata]
        pss_replot(newrange)
    end
    
    #Zoom
    function pss_zoom(zoomsteps::Int64)
        #Calculate scaling (8 steps to scale by factor of 2)
        scale::Float64 = 2.0^(zoomsteps/8)
        
        #Get transformed range
        range::Array{Float64,1} = [plt[:axis]()...]
        centerx::Float64,centery::Float64 = (range[2]+range[1])/2, (range[4]+range[3])/2
        newdx::Float64,newdy::Float64 = scale*(range[2]-range[1])/2, scale*(range[4]-range[3])/2
        newrange::Array{Float64,1} = [centerx-newdx, centerx+newdx,
                                      centery-newdy, centery+newdy]
        pss_replot(newrange)
    end
    
    
    # #########################################################################
    #Event handlers

    #Left clicking selects the mode to be emphasized
    #Right clicking deselects any emphasized mode
    #Middle clicking and dragging moves PSS
    #Scrolling zooms in and out

    #Mouse press-drag-release state
    pss_dragstart::Array{Float64,1} = [0.0,0.0]
    pss_button::Int64 = 1
    pss_pressed::Bool = false

    #Press event
    function onpress(mev::PyObject)
        #Ignore all mouse presses outside axes, or when a button is already pressed
        if !pss_pressed && mev[:inaxes] != nothing
            #Store new press-drag-release state
            pss_pressed = true
            pss_button = mev[:button]
            pss_dragstart = Float64[mev[:xdata],mev[:ydata]]
        end
    end

    #Release event
    function onrelease(mev::PyObject)
        #All button releases other than from the one first pressed are ignored
        if pss_pressed && mev[:button] == pss_button
            #Toggle pressed state
            pss_pressed = false
            
            #Ensure release is in axes
            axes = mev[:inaxes]
            if axes != nothing
                #Choose task depending on button pressed
                if pss_button == 1
                    #Left release
                    #Only clicks without drags are handled
                    if mev[:xdata] == pss_dragstart[1] && mev[:ydata] == pss_dragstart[2]
                        #Select mode with point at click location
                        pss_dehighlight()
                        pss_highlight(mev[:xdata],mev[:ydata])
                    end
                elseif pss_button == 3
                    #Right release
                    #Only clicks without drags are handled
                    if mev[:xdata] == pss_dragstart[1] && mev[:ydata] == pss_dragstart[2]
                        #Deselect mode
                        pss_dehighlight()
                    end
                elseif pss_button == 2
                    #Middle release
                    pss_pan(mev[:xdata]-pss_dragstart[1],mev[:ydata]-pss_dragstart[2])
                end
            end
        end
    end
    
    #Mouse scroll state
    pss_zoomsteps::Int64 = 0
    
    #Scroll event
    function onscroll(mev::PyObject)
        axes = mev[:inaxes]
        if axes != nothing
            #Zoom on scroll
            pss_zoomsteps += mev[:step]
        end
    end
    
    #PSS closed state
    pss_closed::Bool = false
    
    #Close event
    function onclose(cev::PyObject)
        pss_closed = true
    end
    
    # #########################################################################
    #Connect handlers and show plot
    
    #Connect event handlers
    fig[:canvas][:mpl_connect]("button_press_event",onpress)
    fig[:canvas][:mpl_connect]("button_release_event",onrelease)
    fig[:canvas][:mpl_connect]("scroll_event",onscroll)
    fig[:canvas][:mpl_connect]("close_event",onclose)
    
    #Show plot in the background
    plt[:ion]()
    plt[:show]()
    plt[:ioff]()
    
    #Update zoom level at regular intervals
    while !pss_closed
        if pss_zoomsteps != 0
            zoomsteps::Int64 = pss_zoomsteps
            pss_zoomsteps = 0
            pss_zoom(zoomsteps)
        end
        sleep(1)
    end
    
end
