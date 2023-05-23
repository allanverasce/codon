// Decompiled by DJ v3.12.12.101 Copyright 2016 Atanas Neshkov  Date: 27/12/2019 12:13:26
// Home Page:  http://www.neshkov.com/dj.html - Check often for new version!
// Decompiler options: packimports(3) 
// Source File Name:   FeatureTranslator.java

package uk.ac.ebi.kraken.parser.translator;

import com.google.common.base.Strings;
import java.util.ArrayList;
import java.util.List;
import uk.ac.ebi.kraken.interfaces.uniprot.features.*;
import uk.ac.ebi.kraken.model.factories.DefaultFeatureFactory;
import uk.ac.ebi.kraken.parser.FeatureHelper;
import uk.ac.ebi.kraken.parser.Translator;

// Referenced classes of package uk.ac.ebi.kraken.parser.translator:
//            CommentTranslatorHelper

public class FeatureTranslator
    implements Translator
{

    public FeatureTranslator()
    {
    }

    public Feature translate(String annotation)
    {
        int index = annotation.indexOf(" ");
        String type = annotation.substring(0, index);
        annotation = annotation.substring(index + 1).trim();
        Feature feature = createFeature(type);
        if(feature != null)
            translate(annotation, feature);
        return feature;
    }

    public void translate(String annotation, Feature feature)
    {
//    	System.err.println ("<ANO>" + annotation+"</ANO>");
    	return;
    	/*annotation = annotation.replace ("..", " ");
    	System.err.println ("<ANO>" + annotation+"</ANO>");
    	
        if(annotation.startsWith(feature.getType().name()))
            annotation = annotation.substring(feature.getType().name().length() + 1).trim();
        int index = annotation.indexOf("/FTId=");
        String ftid = null;
        if(index != -1)
        {
            ftid = annotation.substring(index + "/FTId=".length());
            if(ftid.endsWith("."))
                ftid = ftid.substring(0, ftid.length() - 1);
            annotation = annotation.substring(0, index).trim();
            if(annotation.endsWith("\n"))
                annotation.substring(0, annotation.length() - 1).trim();
        }
        if(ftid != null && (feature instanceof HasFeatureId))
            ((HasFeatureId)feature).setFeatureId(FACTORY.buildFeatureId(ftid));
        List evidences = new ArrayList();
        String value = CommentTranslatorHelper.stripEvidenceIds(annotation, evidences);
        feature.setEvidenceIds(evidences);
        String locationStart;

    	index = value.indexOf (" ");
    	locationStart = value.substring(0, index);
    	value = value.substring (index + 1);

        String locationEnd;
        index = value.indexOf("\n");
        if(index != -1)
        {
            locationEnd = value.substring (0, index);
            value = value.substring(index + 1);
        } else
        {
        	index = value.indexOf("\n");
            if(index != -1)
            {
                locationEnd = value.substring (0, index);
                value = value.substring(index + 1);
            }
            else {
            	locationEnd = value;
            	value = "";
            }
        }
        updateFeatureLocation(feature, locationStart, locationEnd);
        if(!Strings.isNullOrEmpty(value))
        {
            String text = value.trim();
            if(text.endsWith("."))
                text = text.substring(0, text.length() - 1);
            updateFeatureDescription(feature, text);
        }*/
    }

    public static Feature createFeature(String name)
    {
        FeatureType afeaturetype[] = FeatureType.values();
        int i = afeaturetype.length;
        for(int j = 0; j < i; j++)
        {
            FeatureType featureType = afeaturetype[j];
            if(name.equals(featureType.name()))
                return FACTORY.buildFeature(featureType);
        }

        return null;
    }

    public static void updateFeatureDescription(Feature feature, String text)
    {
        if(feature instanceof ConflictFeature)
            text = FeatureHelper.consumeConflictReport((ConflictFeature)feature, text);
        if((feature instanceof HasFeatureDescription) && !(feature instanceof ConflictFeature) && text != null && !text.isEmpty())
            ((HasFeatureDescription)feature).setFeatureDescription(FACTORY.buildFeatureDescription(text));
        if(feature instanceof MutagenFeature)
            text = FeatureHelper.consumeMutagenReport((MutagenFeature)feature, text);
        if(feature instanceof VariantFeature)
            text = FeatureHelper.consumeVariantReport((VariantFeature)feature, text);
        if(feature instanceof VarSeqFeature)
            text = FeatureHelper.consumeVarSplicFeature((VarSeqFeature)feature, text);
        if(feature instanceof HasAlternativeSequence)
            FeatureHelper.consumeAlternativeSequence((HasAlternativeSequence)feature, text);
        if(feature instanceof CarbohydFeature)
        {
            text = FeatureHelper.consumeCarbohydFeature((CarbohydFeature)feature, text);
            if(text.equals("."))
                text = "";
            if(text.endsWith("."))
                text = text.substring(0, text.length() - 1);
            ((HasFeatureDescription)feature).setFeatureDescription(DefaultFeatureFactory.getInstance().buildFeatureDescription(text.trim()));
        }
    }

    public static void updateFeatureLocation(Feature feature, String locationStart, String locationEnd)
    {
    	System.err.println (locationStart + "  " + locationEnd);
        if(locationStart == null)
            feature.getFeatureLocation().setStartModifier(FeatureLocationModifier.UNKOWN);
        else
        if(locationStart.trim().isEmpty())
        {
            feature.getFeatureLocation().setStartModifier(FeatureLocationModifier.UNKOWN);
        } else
        {
            locationStart = locationStart.trim();
            char c = locationStart.charAt(0);
            if(c == '?')
            {
                if(locationStart.length() > 1)
                {
                    String val = locationStart.substring(1).trim();
                    if(val.isEmpty())
                    {
                        feature.getFeatureLocation().setStartModifier(FeatureLocationModifier.UNKOWN);
                    } else
                    {
                        int value = Integer.parseInt(val);
                        if(value == -1)
                        {
                            feature.getFeatureLocation().setStartModifier(FeatureLocationModifier.UNKOWN);
                        } else
                        {
                            feature.getFeatureLocation().setStart(value);
                            feature.getFeatureLocation().setStartModifier(FeatureLocationModifier.UNSURE);
                        }
                    }
                } else
                {
                    feature.getFeatureLocation().setStartModifier(FeatureLocationModifier.UNKOWN);
                }
            } else
            if(c == '<')
            {
                feature.getFeatureLocation().setStartModifier(FeatureLocationModifier.OUTSIDE_KNOWN_SEQUENCE);
                if(locationStart.length() > 1)
                {
                    String val = locationStart.substring(1);
                    feature.getFeatureLocation().setStart(Integer.parseInt(val.trim()));
                }
            } else
            {
                feature.getFeatureLocation().setStartModifier(FeatureLocationModifier.EXACT);
                feature.getFeatureLocation().setStart(Integer.parseInt(locationStart));
            }
        }
        if(locationEnd == null)
            feature.getFeatureLocation().setEndModifier(FeatureLocationModifier.UNKOWN);
        else
        if(locationEnd.trim().isEmpty())
        {
            feature.getFeatureLocation().setEndModifier(FeatureLocationModifier.UNKOWN);
        } else
        {
            locationEnd = locationEnd.trim();
            char c = locationEnd.charAt(0);
            if(c == '?')
            {
                feature.getFeatureLocation().setEndModifier(FeatureLocationModifier.UNSURE);
                if(locationEnd.length() > 1)
                {
                    String val = locationEnd.substring(1);
                    int value = Integer.parseInt(val);
                    if(value == -1)
                    {
                        feature.getFeatureLocation().setEndModifier(FeatureLocationModifier.UNKOWN);
                    } else
                    {
                        feature.getFeatureLocation().setEnd(value);
                        feature.getFeatureLocation().setEndModifier(FeatureLocationModifier.UNSURE);
                    }
                } else
                {
                    feature.getFeatureLocation().setEndModifier(FeatureLocationModifier.UNKOWN);
                }
            } else
            if(c == '>')
            {
                feature.getFeatureLocation().setEndModifier(FeatureLocationModifier.OUTSIDE_KNOWN_SEQUENCE);
                if(locationEnd.length() > 1)
                {
                    String val = locationEnd.substring(1);
                    feature.getFeatureLocation().setEnd(Integer.parseInt(val.trim()));
                }
            } else
            {
                feature.getFeatureLocation().setEndModifier(FeatureLocationModifier.EXACT);
                feature.getFeatureLocation().setEnd(Integer.parseInt(locationEnd));
            }
        }
    }

//    public volatile Object translate(String s)
//    {
//        return translate(s);
//    }

    public void translate (String s, Object obj)
    {
        translate(s, (Feature)obj);
    }

    private static final String FTID = "/FTId=";
    private static final String SPACE = " ";
    private static final DefaultFeatureFactory FACTORY = DefaultFeatureFactory.getInstance();

}
