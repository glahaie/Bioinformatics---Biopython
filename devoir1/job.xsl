<?xml version="1.0" encoding="ISO-8859-1"?>
<!-- 
job.xsl stylesheet
Transforms a Mobyle xml jobState definition into an HTML form.
This form is either:
- a standalone AJAX-less form that can be called from outside the Mobyle Portal or
- a form (i.e., not a whole html document) that is loaded in the Mobyle Portal
Authors: Hervé Ménager, Bertrand Néron
 -->
<xsl:stylesheet version="1.0"
  xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
  <xsl:output method="html" indent="yes" doctype-system="http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd" doctype-public="-//W3C//DTD XHTML 1.0 Strict//EN"/>

  <!-- if set to 'true', means that we should not generate a complete form, but only one to fill the program channel in the Portal -->
  <xsl:param name="isInPortal" />

  <!-- specifies if the file links are made to the same directory or to a distant server -->
  <xsl:param name="local">yes</xsl:param>

  <xsl:param name="previewDataLimit">1048576</xsl:param>

  <xsl:param name="serviceXmlUri" />
  
  <xsl:variable name="jobStateDoc" select="/"/>
  
  <xsl:include href="mobyle.xsl" />

  <!-- job info in html title -->
  <xsl:template match="//jobState" mode="head">
    <title><xsl:value-of select="concat(name/text(),' job ', id/text())" /></title>    
  </xsl:template>

  <!-- job body -->
  <xsl:template match="//jobState">
    <div class="jobdiv">
      <xsl:attribute name="id">
        <xsl:value-of select="id" />
      </xsl:attribute>
      <xsl:if test="$isInPortal='true'">
        <!-- popups overlay -->
        <xsl:call-template name="popups">
          <xsl:with-param name="popupName" select="id/text()"/>
        </xsl:call-template>
      </xsl:if>      
      <!-- Description -->
      <xsl:if test="$isInPortal!='true'">
        <h2><xsl:value-of select="concat(name/text(),' job ', id/text())" /></h2>
      </xsl:if>
      <input type="hidden">
        <xsl:attribute name="id"><xsl:value-of select="concat(id, '_programName')" /></xsl:attribute>
        <xsl:attribute name="value"><xsl:value-of select="name" /></xsl:attribute>
      </input>
      <input type="hidden">
        <xsl:attribute name="id"><xsl:value-of select="concat(id, '_date')" /></xsl:attribute>
        <xsl:attribute name="value"><xsl:value-of select="date" /></xsl:attribute>
      </input>
      <input type="hidden">
        <xsl:attribute name="id"><xsl:value-of select="concat(id, '_status')" /></xsl:attribute>
        <xsl:attribute name="value"><xsl:value-of select="status/value" /></xsl:attribute>
      </input>
      <input type="hidden">
        <xsl:attribute name="id"><xsl:value-of select="concat(id, '_message')" /></xsl:attribute>
        <xsl:attribute name="value"><xsl:value-of select="status/message" /></xsl:attribute>
      </input>
      
      <div class="header">
        <fieldset class="fieldset">
          <legend>status</legend>
          <table class="resultstable">
            <tr>
              <th>job id</th>
              <td><xsl:value-of select="id" /></td>
              <td>
                <xsl:if test="$isInPortal='true'">
                  <!-- job removal button -->
                  <xsl:call-template name="removeJobButton">
                    <xsl:with-param name="popupName" select="id"/>
                  </xsl:call-template>
                </xsl:if>
              </td>
            </tr>
            <tr>
              <th>program</th>
              <td><xsl:value-of select="name" /></td>
              <td>
                <xsl:if test="$isInPortal='true'">
                  <!-- help popup display button -->
                  <xsl:call-template name="askHelpPopUpButton">
                    <xsl:with-param name="popupName" select="id"/>
                  </xsl:call-template>
                </xsl:if>
              </td>
            </tr>
            <tr>
              <th>submission date</th>
              <td><xsl:value-of select="date" /></td>
              <td>
                <xsl:if test="$isInPortal='true'">
                  <!-- back to the form display button -->
                  <div>
                    <input type="button" value="&lt;&lt; go back to program">
                      <xsl:attribute name="id">
                        <xsl:value-of select="concat('backtoform_',id)" />
                      </xsl:attribute>                    
                    </input>
                  </div>
                </xsl:if>
              </td>
            </tr>
            <tr><th>Mobyle server</th><td><xsl:value-of select="host" /></td></tr>
            <tr>
              <th>status</th>
              <td>
                <xsl:attribute name="class">
                  <xsl:value-of
                    select="concat('jobStatus',status/value)" />
                </xsl:attribute>
                <xsl:value-of select="status/value" />
              </td>
              <td>                
                <xsl:if test="$isInPortal='true'">
                <div><input type="button" class="refresh_link" value="update job status" /></div>
                </xsl:if>      
              </td>
            </tr>
            <xsl:if test="status/message">
              <tr>
                <th>message</th>
                <td class="errormsg">
                  <xsl:value-of select="status/message" />
                </td>
              </tr>
            </xsl:if>
          </table>
        </fieldset>
      </div>
      <div class="result">
        <xsl:if test="data/output">
        <fieldset class="fieldset">
          <legend>job results</legend>
		  <xsl:variable name="documentAvailable" select="function-available('document')"/>	
          <xsl:choose>
            <xsl:when test="$documentAvailable='true' and $serviceXmlUri">
			  <xsl:variable name="serviceDoc" select="document($serviceXmlUri)"/>
              <xsl:apply-templates select="$serviceDoc//program/parameters"/>
              <xsl:if test="(//output/parameter/name/text()='stdout') and not($serviceDoc//program/parameters/@isstdout)">
                <xsl:apply-templates select="data/output[parameter/name/text()='stdout']"/>
              </xsl:if>
              <xsl:if test="//output/parameter/name/text()='stderr'">
                <xsl:apply-templates select="data/output[parameter/name/text()='stderr']"/>              
              </xsl:if>
            </xsl:when>
            <xsl:otherwise>
              <xsl:apply-templates select="data/output"/>            
            </xsl:otherwise>
          </xsl:choose>
        </fieldset>
        </xsl:if>
      </div>
      <div class="param">
        <fieldset class="fieldset">
          <legend>parameters</legend>
          <xsl:apply-templates select="data/input"/>
          <xsl:apply-templates select="commandLine"/>
          <xsl:apply-templates select="paramFiles"/>
        </fieldset>
      </div>
      <xsl:if test="$isInPortal='true' and (status/value='finished')">
        <fieldset class="fieldset zipfile">
          <legend>job archive</legend>
          <a target="_blank">
           <xsl:attribute name="href">
            <xsl:value-of select="concat(id/text(), '/', name/text(),'_',substring-after(id/text(),concat(name/text(),'/')),'.zip')" />
           </xsl:attribute>
           download this job as an archive
          </a>
        </fieldset>
      </xsl:if>
    </div>  
  </xsl:template>

  <!-- template for each parameter node -->
  <xsl:template match="parameter">
    <xsl:if test="(boolean(@isout) or boolean(@isstdout)) and not(boolean(@ishidden))">
      <div class="parameter">
        <xsl:apply-templates select="$jobStateDoc//output[parameter/name/text()=current()/name/text()]"/>
      </div>
    </xsl:if>
  </xsl:template>
  
  <xsl:template match="paragraph">
    <xsl:if test=".//parameter[(boolean(@isout) or boolean(@isstdout)) and ./name/text()=$jobStateDoc//output/parameter/name/text()]">
      <fieldset class="fieldset">
        <legend>
          <xsl:value-of select="prompt" />
          <xsl:call-template name="commentToggle">
            <xsl:with-param name="commentNode" select="comment"/>
          </xsl:call-template>
        </legend>
        <xsl:call-template name="commentText">
          <xsl:with-param name="commentNode" select="comment"/>
        </xsl:call-template>
        <xsl:choose>
          <xsl:when test="not(layout)">
            <xsl:apply-templates select="parameters" />
          </xsl:when>
          <xsl:otherwise>
            <xsl:apply-templates select="layout" />            
          </xsl:otherwise>
        </xsl:choose>
      </fieldset>
    </xsl:if>
  </xsl:template>

  <!-- job output, no order -->
  <xsl:template match="output">
    <xsl:if test="string-length(file/text())>0 or string-length(value/text())>0">
    <!-- output title -->
    <fieldset class="fieldset">
      <legend>
        <xsl:choose>
          <xsl:when test="parameter/prompt">
            <xsl:value-of select="parameter/prompt" /><!-- prompt -->
          </xsl:when>
          <xsl:otherwise>
            <xsl:value-of select="parameter/name" /><!-- parameter name -->
          </xsl:otherwise>
        </xsl:choose>
      (<xsl:value-of select="normalize-space(concat(parameter/type/biotype, ' ', parameter/type/datatype/class))" />)
      <xsl:call-template name="commentToggle">
        <xsl:with-param name="commentNode" select="parameter/comment"/>
        <xsl:with-param name="exampleNode" select="parameter/example"/>
      </xsl:call-template>        
      </legend>
      <xsl:call-template name="commentText">
        <xsl:with-param name="commentNode" select="parameter/comment"/>
        <xsl:with-param name="exampleNode" select="parameter/example"/>
      </xsl:call-template>
      <!-- data -->
      <xsl:apply-templates select="file" mode="preview"/>
    </fieldset>
    </xsl:if>
  </xsl:template>

  <!-- job input -->
  <xsl:template match="input">
    <fieldset class="fieldset">
      <legend>
        <xsl:choose>
          <xsl:when test="parameter/prompt">
            <xsl:value-of select="parameter/prompt" /><!-- prompt -->
          </xsl:when>
          <xsl:otherwise>
            <xsl:value-of select="parameter/name" /><!-- parameter name -->
          </xsl:otherwise>
        </xsl:choose>
      (<xsl:value-of select="normalize-space(concat(parameter/type/biotype, ' ', parameter/type/datatype/class))" />)
      <xsl:call-template name="commentToggle">
        <xsl:with-param name="commentNode" select="parameter/comment"/>
        <xsl:with-param name="exampleNode" select="parameter/example"/>
      </xsl:call-template>        
      </legend>
      <xsl:call-template name="commentText">
        <xsl:with-param name="commentNode" select="parameter/comment"/>
        <xsl:with-param name="exampleNode" select="parameter/example"/>
      </xsl:call-template>
      <xsl:if test="file">
        <table border="1">
          <th>file</th><th>format detection program</th><th>re-formatted file</th>
          <tr>
            <td>
              <xsl:apply-templates select="file" mode="link"/>
            </td>
            <td>
              <xsl:apply-templates select="fmtProgram"/>
            </td>
            <td>
              <xsl:apply-templates select="formattedFile" mode="link"/>
            </td>
          </tr>
        </table>
      </xsl:if>
      <xsl:if test="value">
        <xsl:apply-templates select="value" />
      </xsl:if>
    </fieldset>
  </xsl:template>

  <!--  job file (input or output) displayed as a link -->
  <xsl:template match="file | formattedFile" mode="link">
    <span>
      <a target="_blank">
        <xsl:variable name='fileLink' select="concat(/jobState/id,'/',.)"/>
        <xsl:attribute name="href"><xsl:value-of select="$fileLink" /></xsl:attribute>
        <xsl:value-of select="text()" />
      </a>
      <xsl:if test="@fmt">
        (<xsl:value-of select="@fmt" /> format)
      </xsl:if>
    </span>
  </xsl:template>

  <!--  job input value -->
  <xsl:template match="value">
    <span>
      Value: <xsl:value-of select="text()" />
    </span>
  </xsl:template>

  <!--  job file (input or output) displayed as a preview -->
  <xsl:template match="file" mode="preview">
    <div class="results_file">
      <xsl:variable name='resultId' select="concat(/jobState/id,'/',.)"/>
      <input type="hidden" class="parameterName">
        <xsl:attribute name="name"><xsl:value-of select="concat($resultId, '_parameterName')" /></xsl:attribute>
        <xsl:attribute name="value"><xsl:value-of select="../parameter/name" /></xsl:attribute>
      </input>
      <input type="hidden" class="fileName">
        <xsl:attribute name="name"><xsl:value-of select="concat($resultId, '_fileName')" /></xsl:attribute>
        <xsl:attribute name="value"><xsl:value-of select="./text()" /></xsl:attribute>
      </input>
      <input type="hidden" class="resultId">
        <xsl:attribute name="name"><xsl:value-of select="concat($resultId, '_resultId')" /></xsl:attribute>
        <xsl:attribute name="value"><xsl:value-of select="$resultId" /></xsl:attribute>
      </input>
      <input type="hidden" class="parameterName">
        <xsl:attribute name="name"><xsl:value-of select="concat($resultId, '_parameterName')" /></xsl:attribute>
        <xsl:attribute name="value"><xsl:value-of select="../parameter/name" /></xsl:attribute>
      </input>
      <input type="hidden" class="parameterDataType">
        <xsl:attribute name="name"><xsl:value-of select="concat($resultId, '_parameterDataType')" /></xsl:attribute>
        <xsl:attribute name="value"><xsl:value-of select="../parameter/type/datatype/class" /></xsl:attribute>
      </input>
      <xsl:for-each select="../parameter/type/biotype">
        <input type="hidden" class="parameterBioType">
          <xsl:attribute name="name"><xsl:value-of select="concat($resultId, '_parameterBioType_', position())" /></xsl:attribute>
          <xsl:attribute name="value"><xsl:value-of select="." /></xsl:attribute>
        </input>
      </xsl:for-each>
      <input type="hidden" class="parameterCard">
        <xsl:attribute name="name"><xsl:value-of select="concat($resultId, '_parameterCard')" /></xsl:attribute>
        <xsl:attribute name="value"><xsl:value-of select="../parameter/type/card" /></xsl:attribute>
      </input>
      <div>
        <xsl:value-of select="." />
		<a target="_blank" class="saveFileLink" title="save this file">
	        <xsl:attribute name="href"><xsl:value-of select="$resultId" />?save</xsl:attribute>
			save
		</a>
      </div><!-- result file name -->
      <xsl:variable name='fileLink' select="concat(/jobState/id,'/',.)"/>
	  <xsl:choose>
	  	<xsl:when test="./@size&lt;=$previewDataLimit">
	      <a target="_blank" title="display in a new window">
	        <xsl:attribute name="href"><xsl:value-of select="$resultId" /></xsl:attribute>
	        <object>
	          <xsl:attribute name="data">
	            <xsl:value-of select="$fileLink" />
	          </xsl:attribute>
	          <xsl:attribute name="type">
	            <xsl:variable name="extension" select="substring-after(.,'.')"/>
	            <xsl:choose>
	              <xsl:when test="($extension = 'out') or ($extension = 'err')">text/plain</xsl:when>
	              <xsl:when test="$extension = 'html'">text/html</xsl:when>
	              <xsl:when test="($extension = 'ps') or ($extension = 'eps')">application/postscript</xsl:when>
	              <xsl:when test="($extension = 'pdf')">application/postscript</xsl:when>
	              <xsl:when test="$extension = 'gif'">image/gif</xsl:when>
	              <xsl:when test="$extension = 'png'">image/png</xsl:when>
	              <xsl:otherwise>text/plain</xsl:otherwise>
	            </xsl:choose>
	          </xsl:attribute>
	        </object>
	      </a>
	    </xsl:when>
		<xsl:otherwise>
			<p class="commentText">The file is too big to be safely displayed here.</p>
		</xsl:otherwise>
	  </xsl:choose>
      <!-- result box -->
      <xsl:apply-templates select="." mode="resultbox"/>
    </div>
  </xsl:template>  

  <xsl:template match="file" mode="resultbox">
    <xsl:variable name='resultId' select="concat(/jobState/id,'/',.)"/>
    <div class="resultbox">
        <span>
          <form target="_blank">
            <xsl:attribute name="action"><xsl:value-of select="$resultId" /></xsl:attribute>
            <input type="submit" value="full screen view" title="display in a new window" />
          </form>
        </span>
        <xsl:if test="$isInPortal='true'">
          <span>        
            result name: 
            <input size="8" type="text">
              <xsl:attribute name="name"><xsl:value-of select="$resultId" /></xsl:attribute>
              <xsl:attribute name="id"><xsl:value-of select="$resultId" /></xsl:attribute>
              <xsl:attribute name="value"><xsl:value-of select="." /></xsl:attribute>
            </input>
          </span>
          <span>        
            <input type="button" class="saveresultlink" title="bookmark this result" value="bookmark"/>
            <span class="pipebox">
               or 
              <select class="pipeselect">
              </select>
              <input type="button" class="piperesultlink" href="#" title="use this file as input data in a new job" value="further analysis" />
            </span>
          </span>
        </xsl:if>
      </div>  
  </xsl:template>

  <!--  job command line -->
  <xsl:template match="commandLine">
    <fieldset class="fieldset">
      <legend>Command line</legend>
        <xsl:value-of select="text()" />
    </fieldset>
  </xsl:template>

  <!--  parameters file -->
  <xsl:template match="paramFiles">
    <fieldset class="fieldset">
      <legend>Parameter file(s)</legend>
	  <ul class="paramFiles">
      <xsl:for-each select="file">
      	<li>
        <a target="_blank">
          <xsl:attribute name="href">
            <xsl:value-of select="concat(/jobState/id,'/',text())" />
          </xsl:attribute>
          <xsl:value-of select="text()" />
        </a>
		</li>
      </xsl:for-each>
	  </ul>
    </fieldset>
  </xsl:template>

</xsl:stylesheet>
